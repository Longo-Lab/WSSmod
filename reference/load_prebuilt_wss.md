# Load a Prebuilt WSS Module/Weighting Table

Reads one of the module weighting tables shipped with WSSmod and applies
the same deduplication and filtering used to build the reference sets:
for proteins assigned to multiple modules, only the assignment with the
highest `mean_alpha_scaled` is kept, and modules with two or fewer
proteins are dropped.

## Usage

``` r
load_prebuilt_wss(prebuilt = "core_AD_plasma_biomarkers")
```

## Arguments

- prebuilt:

  Name of a prebuilt reference set. See
  [`list_prebuilt_wss()`](https://Longo-Lab.github.io/WSSmod/reference/list_prebuilt_wss.md)
  for available options.

## Value

A data.table with columns `module`, `symbol`, `mean_beta`, and
`mean_alpha_scaled`, suitable for use as the `result` argument of
[`calculate_WSS()`](https://Longo-Lab.github.io/WSSmod/reference/calculate_WSS.md).

## See also

[`calculate_WSS()`](https://Longo-Lab.github.io/WSSmod/reference/calculate_WSS.md),
[`list_prebuilt_wss()`](https://Longo-Lab.github.io/WSSmod/reference/list_prebuilt_wss.md),
[`wss_prebuilt_terms()`](https://Longo-Lab.github.io/WSSmod/reference/wss_prebuilt_terms.md)
