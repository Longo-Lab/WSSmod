# Predict WSS Module Scores from Core Biomarkers

Predicts module scores for a prebuilt reference set's `joinet` model (a
two-layer elastic net stack) without needing to measure the full
proteomic panel. For `"core_AD_plasma_biomarkers"`, this predicts all 75
module scores from Age, Gender, and the 8 core plasma biomarkers used in
the original analysis.

## Usage

``` r
predict_WSS(newdata, prebuilt = "core_AD_plasma_biomarkers", type = "response")
```

## Arguments

- newdata:

  A data.frame or matrix with one row per sample and (at least) the
  columns required by the model; see `x_cols` in
  [`load_prebuilt_wss_model()`](https://Longo-Lab.github.io/WSSmod/reference/load_prebuilt_wss_model.md)
  for the exact set and order. For `"core_AD_plasma_biomarkers"` this is
  `Age` (raw years), `Gender` (numeric, `1` = male, `0` = not male), and
  8 biomarker columns suffixed `_norm` (`PlasmaPTau181_norm`,
  `PlasmaAB142P_norm`, `PlasmaAB140P_norm`, `PlasmaABRatio_norm`,
  `PlasmapTau217_norm`, `PlasmapTau217_AB42Ratio_norm`,
  `PlasmaGFAP_norm`, `PlasmaNfL_norm`).

  **The `_norm` biomarker columns must already be rank-based
  inverse-normal transformed (e.g. via `RNOmni::RankNorm()`) relative to
  your own reference cohort before calling this function.** This
  transform is inherently relative to the distribution it's computed
  against, so it cannot be done correctly for a single new sample in
  isolation, and this function does not attempt it for you.

- prebuilt:

  Name of a prebuilt reference set. See
  [`list_prebuilt_wss()`](https://Longo-Lab.github.io/WSSmod/reference/list_prebuilt_wss.md)
  for available options.

- type:

  Prediction type passed to `predict.joinet()`: `"response"` (the
  default) or `"link"`.

## Value

A list with components `base` (first-layer-only predictions) and `meta`
(final stacked predictions), each a matrix with one row per sample (row
names taken from `newdata`, if present) and one column per module (named
by the model's `outcomes`).

## See also

[`load_prebuilt_wss_model()`](https://Longo-Lab.github.io/WSSmod/reference/load_prebuilt_wss_model.md),
[`calculate_WSS()`](https://Longo-Lab.github.io/WSSmod/reference/calculate_WSS.md)

## Examples

``` r
if (FALSE) { # \dontrun{
newdata <- data.frame(
  Age = 72,
  Gender = 1,
  PlasmaPTau181_norm = 0,
  PlasmaAB142P_norm = 0,
  PlasmaAB140P_norm = 0,
  PlasmaABRatio_norm = 0,
  PlasmapTau217_norm = 0,
  PlasmapTau217_AB42Ratio_norm = 0,
  PlasmaGFAP_norm = 0,
  PlasmaNfL_norm = 0
)
pred <- predict_WSS(newdata)
pred$meta
} # }
```
