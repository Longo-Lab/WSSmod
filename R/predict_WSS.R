#' Predict WSS Module Scores from Core Biomarkers
#'
#' Predicts module scores for a prebuilt reference set's `joinet` model (a
#' two-layer elastic net stack) without needing to measure the full
#' proteomic panel. For `"core_AD_plasma_biomarkers"`, this predicts all 75
#' module scores from Age, Gender, and the 8 core plasma biomarkers used in
#' the original analysis.
#'
#' @param newdata A data.frame or matrix with one row per sample and (at
#'   least) the columns required by the model; see `x_cols` in
#'   [load_prebuilt_wss_model()] for the exact set and order. For
#'   `"core_AD_plasma_biomarkers"` this is `Age` (raw years), `Gender`
#'   (numeric, `1` = male, `0` = not male), and 8 biomarker columns suffixed
#'   `_norm` (`PlasmaPTau181_norm`, `PlasmaAB142P_norm`, `PlasmaAB140P_norm`,
#'   `PlasmaABRatio_norm`, `PlasmapTau217_norm`,
#'   `PlasmapTau217_AB42Ratio_norm`, `PlasmaGFAP_norm`, `PlasmaNfL_norm`).
#'
#'   **The `_norm` biomarker columns must already be rank-based
#'   inverse-normal transformed (e.g. via [RNOmni::RankNorm()]) relative to
#'   your own reference cohort before calling this function.** This
#'   transform is inherently relative to the distribution it's computed
#'   against, so it cannot be done correctly for a single new sample in
#'   isolation, and this function does not attempt it for you.
#' @param prebuilt Name of a prebuilt reference set. See
#'   [list_prebuilt_wss()] for available options.
#' @param type Prediction type passed to `predict.joinet()`: `"response"`
#'   (the default) or `"link"`.
#'
#' @return A list with components `base` (first-layer-only predictions) and
#'   `meta` (final stacked predictions), each a matrix with one row per
#'   sample (row names taken from `newdata`, if present) and one column per
#'   module (named by the model's `outcomes`).
#'
#' @seealso [load_prebuilt_wss_model()], [calculate_WSS()]
#'
#' @examples
#' \dontrun{
#' newdata <- data.frame(
#'   Age = 72,
#'   Gender = 1,
#'   PlasmaPTau181_norm = 0,
#'   PlasmaAB142P_norm = 0,
#'   PlasmaAB140P_norm = 0,
#'   PlasmaABRatio_norm = 0,
#'   PlasmapTau217_norm = 0,
#'   PlasmapTau217_AB42Ratio_norm = 0,
#'   PlasmaGFAP_norm = 0,
#'   PlasmaNfL_norm = 0
#' )
#' pred <- predict_WSS(newdata)
#' pred$meta
#' }
#'
#' @export
predict_WSS <- function(newdata, prebuilt = "core_AD_plasma_biomarkers", type = "response") {
  model_info <- load_prebuilt_wss_model(prebuilt)
  if (is.null(model_info)) {
    stop(sprintf(
      "Prebuilt reference set '%s' has no associated prediction model", prebuilt
    ), call. = FALSE)
  }

  missing_cols <- setdiff(model_info$x_cols, colnames(newdata))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "newdata is missing required column(s): %s",
      paste(missing_cols, collapse = ", ")
    ), call. = FALSE)
  }

  newdata <- newdata[, model_info$x_cols, drop = FALSE]
  non_numeric <- names(newdata)[!vapply(newdata, is.numeric, logical(1))]
  if (length(non_numeric) > 0) {
    stop(sprintf(
      "newdata column(s) must be numeric: %s",
      paste(non_numeric, collapse = ", ")
    ), call. = FALSE)
  }

  for (pkg in c("glmnet", "joinet")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf(
        "Package '%s' is required for predict_WSS() but is not installed", pkg
      ), call. = FALSE)
    }
    loadNamespace(pkg)
  }

  newx <- as.matrix(newdata)
  pred <- stats::predict(model_info$model, newx = newx, type = type)

  colnames(pred$base) <- colnames(pred$meta) <- model_info$outcomes
  rownames(pred$base) <- rownames(pred$meta) <- rownames(newx)

  pred
}
