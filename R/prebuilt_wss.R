# Registry of prebuilt WSS reference sets shipped in inst/extdata.
# Each entry names the scores file (module/protein weighting table, required)
# and an optional terms file (human-readable module labels).
.wss_prebuilt_registry <- list(
  core_AD_plasma_biomarkers = list(
    scores = "core_AD_plasma_biomarkers.STRING_merged_scaled_scores.csv",
    terms = "core_AD_plasma_biomarkers.STRING_INT_WSS_terms_Zsum.csv"
  )
)

.wss_prebuilt_entry <- function(prebuilt) {
  entry <- .wss_prebuilt_registry[[prebuilt]]
  if (is.null(entry)) {
    stop(sprintf(
      "Unknown prebuilt WSS reference set '%s'. Available sets: %s",
      prebuilt,
      paste(names(.wss_prebuilt_registry), collapse = ", ")
    ), call. = FALSE)
  }
  entry
}

.wss_prebuilt_path <- function(file) {
  path <- system.file("extdata", file, package = "WSSmod")
  if (!nzchar(path)) {
    stop(sprintf("Could not find packaged reference file '%s'", file), call. = FALSE)
  }
  path
}

#' Load a Prebuilt WSS Module/Weighting Table
#'
#' Reads one of the module weighting tables shipped with WSSmod. The raw
#' table may assign the same protein to more than one module (e.g. once per
#' biomarker it was associated with); [calculate_WSS()] deduplicates and
#' filters this table before use via `dedup_by` and `min_module_size`.
#'
#' @param prebuilt Name of a prebuilt reference set. See
#'   [list_prebuilt_wss()] for available options.
#'
#' @return A data.table with (at least) columns `module`, `symbol`,
#'   `mean_beta`, and `mean_alpha_scaled`, suitable for use as the `result`
#'   argument of [calculate_WSS()].
#'
#' @seealso [calculate_WSS()], [list_prebuilt_wss()], [wss_prebuilt_terms()]
#'
#' @importFrom data.table fread
#' @export
load_prebuilt_wss <- function(prebuilt = "core_AD_plasma_biomarkers") {
  entry <- .wss_prebuilt_entry(prebuilt)

  result <- fread(.wss_prebuilt_path(entry$scores))
  result[, module := paste(Biomarker, cluster, sep = ".")]

  result[]
}

#' List Available Prebuilt WSS Reference Sets
#'
#' @return A character vector of prebuilt reference set names usable as the
#'   `prebuilt` argument of [calculate_WSS()] and [load_prebuilt_wss()].
#'
#' @export
list_prebuilt_wss <- function() {
  names(.wss_prebuilt_registry)
}

#' Get Module Term Labels for a Prebuilt WSS Reference Set
#'
#' Some prebuilt reference sets ship a companion table of human-readable
#' term labels for each module, summarized from the underlying pathway
#' enrichment used to build the module.
#'
#' @param prebuilt Name of a prebuilt reference set. See
#'   [list_prebuilt_wss()] for available options.
#'
#' @return A data.table of module term labels, or `NULL` if the given
#'   reference set has no associated terms file.
#'
#' @importFrom data.table fread
#' @export
wss_prebuilt_terms <- function(prebuilt = "core_AD_plasma_biomarkers") {
  entry <- .wss_prebuilt_entry(prebuilt)

  if (is.null(entry$terms)) {
    return(NULL)
  }

  fread(.wss_prebuilt_path(entry$terms))
}
