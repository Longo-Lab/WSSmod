#' Provenance of the core_AD_plasma_biomarkers files in inst/extdata
#'
#' Records how each file of the "core_AD_plasma_biomarkers" prebuilt reference
#' set was made. The upstream analysis lives in the private Longo-Lab
#' repository Plasma_Ptau_Proteomics (Box: plasma_ptau_proteomics), and the
#' participant-level inputs are not distributed, so this script can only be run
#' by the Longo Lab with Box mounted. It is excluded from the built package
#' (see .Rbuildignore).
#'
#' Files written to `out_dir` (default `inst/extdata`):
#'
#' * `STRING_merged_scaled_scores.csv`: module/protein weighting table, copied
#'   unchanged from the paper pipeline.
#' * `STRING_INT_WSS_terms_Zsum.csv`: module term labels and `mean.Zsum`
#'   (Zsummary preservation across 12 independent GNPC cohorts), copied
#'   unchanged from the paper pipeline.
#' * `prediction_model.rds`: the `Biomarkers_Age_Sex` joinet model, extracted
#'   from the STRING_WSS_v2 final models (Plasma_Ptau_Proteomics#76).
#' * `baseline_biomarker_reference.csv`: anonymized baseline biomarker
#'   reference. Real values from the 389 baseline samples in the joint
#'   complete-case set the model was fit on, with each biomarker independently
#'   shuffled (seed 42) so no values can be linked across biomarkers.
#'   Only the eight biomarker columns are written. The outcomes, Age, Gender,
#'   genetic PCs, APOE and HybControlNormScale are loaded solely to define the
#'   joint complete-case filter and are not shipped.
#'
#' data.table is used for the reference because its merge order fixes the row
#' order that the seeded shuffle acts on; changing it changes the output.
#'
#' Butler RR et al. (2026) <doi:10.64898/2026.08.14.26360394>
#'
#' @author Robert B Butler III

library(data.table)

box_dir <- "/mnt/box/Robert Butler's Files/plasma_ptau_proteomics"
if (!exists("out_dir")) out_dir <- "inst/extdata"
prefix <- "core_AD_plasma_biomarkers"

#' Path of an output file in the reference set
#'
#' @param file File name without the reference set prefix.
#'
#' @return A path inside `out_dir`.
out_path <- function(file) {
  file.path(out_dir, paste(prefix, file, sep = "."))
}

biomarkers <- c(
  "PlasmaPTau181", "PlasmaAB142P", "PlasmaAB140P", "PlasmaABRatio",
  "PlasmapTau217", "PlasmapTau217_AB42Ratio", "PlasmaGFAP", "PlasmaNfL"
)

# --- Weighting table and term labels: unchanged copies ----------------------
file.copy(
  file.path(box_dir, "module_scores/STRING_merged_scaled_scores.csv"),
  out_path("STRING_merged_scaled_scores.csv"),
  overwrite = TRUE
)
file.copy(
  file.path(box_dir, "module_scores/STRING.INT.WSS.terms.Zsum.csv"),
  out_path("STRING_INT_WSS_terms_Zsum.csv"),
  overwrite = TRUE
)

# --- Prediction model: the Biomarkers_Age_Sex variant ------------------------
final_models <- readRDS(file.path(
  box_dir, "prediction_models/STRING_WSS_v2/model_info/joinet/final_models_IN.rds"
))
model_info <- final_models[[2]]
stopifnot(model_info$variant == "Biomarkers_Age_Sex")
saveRDS(model_info, out_path("prediction_model.rds"))

# --- Anonymized baseline biomarker reference --------------------------------
set.seed(42)

baseline <- fread(file.path(box_dir, "dataProt_SS-205063.merged.QC.csv"))
pcs <- fread(
  file.path(box_dir, "genomic/mypcair.PCs.csv"),
  colClasses = c(PIDN = "character")
)
baseline <- merge(baseline, pcs, by = "PIDN")

baseline_wss <- fread(
  file.path(
    box_dir,
    "module_scores/dataProt_SS-205063.STRING.INT.baseline.WSS.csv"
  ),
  colClasses = c(SampleID = "character")
)
baseline <- merge(baseline, baseline_wss, by = "SampleID")

baseline[, Gender := as.numeric(Gender == "M")]
baseline[, APOE.e34 := as.numeric(APOE == "e34")]
baseline[, APOE.e44 := as.numeric(APOE == "e44")]

# Joint complete cases across every outcome, covariate and biomarker column
# used to fit the model: this is the pipeline's frozen reference (389 rows),
# not a per-biomarker NA drop
outcomes <- setdiff(names(baseline_wss), "SampleID")
covariates <- c(
  "Age", "Gender", "HybControlNormScale",
  grep("^PC\\d{1,2}$", names(baseline), value = TRUE),
  "APOE.e34", "APOE.e44"
)
model_cols <- c("SampleID", outcomes, covariates, biomarkers)
reference <- baseline[complete.cases(baseline[, ..model_cols]), ..model_cols]
stopifnot(nrow(reference) == 389)

# Independent per-biomarker shuffle, in the order of `biomarkers`
reference_long <- rbindlist(lapply(biomarkers, function(b) {
  data.table(biomarker = b, value = as.numeric(sample(reference[[b]])))
}))
stopifnot(nrow(reference_long) == 389 * length(biomarkers))

write.csv(
  reference_long,
  out_path("baseline_biomarker_reference.csv"),
  row.names = FALSE
)
