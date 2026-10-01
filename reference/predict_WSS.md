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
  the 8 core biomarkers (`PlasmaPTau181`, `PlasmaAB142P`,
  `PlasmaAB140P`, `PlasmaABRatio`, `PlasmapTau217`,
  `PlasmapTau217_AB42Ratio`, `PlasmaGFAP`, `PlasmaNfL`).

  **The biomarker columns must already be rank-based inverse-normal
  transformed before calling this function** – see
  [`normalize_wss_biomarkers()`](https://Longo-Lab.github.io/WSSmod/reference/normalize_wss_biomarkers.md)
  to build them from raw biomarker values, including for a single new
  patient.

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

[`normalize_wss_biomarkers()`](https://Longo-Lab.github.io/WSSmod/reference/normalize_wss_biomarkers.md),
[`load_prebuilt_wss_model()`](https://Longo-Lab.github.io/WSSmod/reference/load_prebuilt_wss_model.md),
[`calculate_WSS()`](https://Longo-Lab.github.io/WSSmod/reference/calculate_WSS.md)

## Examples

``` r
if (requireNamespace("glmnet", quietly = TRUE) &&
    requireNamespace("joinet", quietly = TRUE)) {
  # Simulated raw biomarkers; works for a single patient too
  raw_biomarkers <- wss_example_biomarkers[, -(1:2)]
  normalized <- normalize_wss_biomarkers(raw_biomarkers)

  newdata <- cbind(
    Age = wss_example_biomarkers$Age,
    Gender = wss_example_biomarkers$Gender,
    normalized
  )
  pred <- predict_WSS(newdata)
  pred$meta[1:3, 1:4]
}
#>           Merged.M17  Merged.M7  Merged.M4 PlasmaPTau181.M26
#> patient01  0.0699344 -0.1617110  0.1312342        -0.2210809
#> patient02  0.1919225  0.0521312 -0.6218288         0.2495096
#> patient03 -0.1928722  0.2029449  0.4536169         0.2990079
```
