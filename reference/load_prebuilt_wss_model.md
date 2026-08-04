# Load a Prebuilt WSS Prediction Model

Some prebuilt reference sets ship a fitted `joinet` model (a two-layer
elastic net) that predicts all of that reference set's module scores
from a small set of covariates, without needing to measure the full
proteomic panel. See
[`predict_WSS()`](https://Longo-Lab.github.io/WSSmod/reference/predict_WSS.md)
to predict from this model directly.

## Usage

``` r
load_prebuilt_wss_model(prebuilt = "core_AD_plasma_biomarkers")
```

## Arguments

- prebuilt:

  Name of a prebuilt reference set. See
  [`list_prebuilt_wss()`](https://Longo-Lab.github.io/WSSmod/reference/list_prebuilt_wss.md)
  for available options.

## Value

A list with components `model` (the fitted `joinet` object), `outcomes`
(the module names it predicts), `x_cols` (the required predictor
columns, in order), and `variant`; or `NULL` if the given reference set
has no associated prediction model.

## See also

[`predict_WSS()`](https://Longo-Lab.github.io/WSSmod/reference/predict_WSS.md),
[`list_prebuilt_wss()`](https://Longo-Lab.github.io/WSSmod/reference/list_prebuilt_wss.md)
