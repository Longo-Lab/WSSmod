# Simulates the example datasets shipped in data/. Both are fully synthetic
# (random draws with a fixed seed); they contain no participant data. Run from
# the package root with devtools::load_all() so the internal helpers are found.

devtools::load_all()
set.seed(20261001)

# --- wss_example_expression ---------------------------------------------------
# Samples x proteins matrix whose columns are the proteins that survive
# calculate_WSS()'s default dedup/size filter on the bundled reference set, so
# calculate_WSS(wss_example_expression) runs without "no proteins" warnings.
# Each module gets a latent per-sample factor; proteins load on it in the
# direction of their mean_beta, plus independent noise.
result <- WSSmod:::.filter_wss_result(load_prebuilt_wss())
n_samples <- 50
sample_ids <- sprintf("sample%02d", seq_len(n_samples))

modules <- unique(result$module)
latent <- matrix(
  rnorm(n_samples * length(modules)),
  nrow = n_samples,
  dimnames = list(sample_ids, modules)
)

expr <- vapply(seq_len(nrow(result)), function(i) {
  loading <- 0.5 * sign(result$mean_beta[i])
  loading * latent[, result$module[i]] + rnorm(n_samples)
}, numeric(n_samples))
dimnames(expr) <- list(sample_ids, result$symbol)

# Rounded to keep the shipped file small
wss_example_expression <- round(expr, 3)

# --- wss_example_biomarkers ---------------------------------------------------
# Age, Gender (1 = male) and the 8 raw plasma biomarkers, drawn from
# log-normal distributions at plausible scales.
n_patients <- 20
medians <- c(
  PlasmaPTau181 = 1.7, PlasmaAB142P = 28, PlasmaAB140P = 305, PlasmaABRatio = 0.092,
  PlasmapTau217 = 0.2, PlasmapTau217_AB42Ratio = 0.007, PlasmaGFAP = 62, PlasmaNfL = 21
)
sdlog <- c(0.4, 0.25, 0.2, 0.2, 0.8, 0.9, 0.5, 0.5)

biomarkers <- mapply(
  function(med, s) round(rlnorm(n_patients, meanlog = log(med), sdlog = s), 4),
  medians, sdlog
)

wss_example_biomarkers <- data.frame(
  Age = round(runif(n_patients, 55, 90)),
  Gender = rbinom(n_patients, 1, 0.5),
  biomarkers,
  row.names = sprintf("patient%02d", seq_len(n_patients))
)

usethis::use_data(
  wss_example_expression, wss_example_biomarkers,
  compress = "xz", overwrite = TRUE
)
