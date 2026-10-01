# Load a Prebuilt WSS Module/Weighting Table

Reads one of the module weighting tables shipped with WSSmod. The raw
table may assign the same protein to more than one module (e.g. once per
biomarker it was associated with);
[`calculate_WSS()`](https://Longo-Lab.github.io/WSSmod/reference/calculate_WSS.md)
deduplicates and filters this table before use via `dedup_by` and
`min_module_size`.

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

A data.table with (at least) columns `module`, `symbol`, `mean_beta`,
and `mean_alpha_scaled`, suitable for use as the `result` argument of
[`calculate_WSS()`](https://Longo-Lab.github.io/WSSmod/reference/calculate_WSS.md).

## See also

[`calculate_WSS()`](https://Longo-Lab.github.io/WSSmod/reference/calculate_WSS.md),
[`list_prebuilt_wss()`](https://Longo-Lab.github.io/WSSmod/reference/list_prebuilt_wss.md),
[`wss_prebuilt_terms()`](https://Longo-Lab.github.io/WSSmod/reference/wss_prebuilt_terms.md)

## Examples

``` r
ref <- load_prebuilt_wss()
head(ref[, c("module", "symbol", "mean_beta", "mean_alpha_scaled")])
#>        module   symbol  mean_beta mean_alpha_scaled
#>        <char>   <char>      <num>             <num>
#> 1: Merged.M17   ZWILCH -0.4732375         1.0000000
#> 2:  Merged.M7    HYOU1  0.1603553         0.9961153
#> 3:  Merged.M4    EFNB2  0.4591349         0.8693081
#> 4:  Merged.M7  DNAJB11  0.2211760         0.8307261
#> 5:  Merged.M7 HSP90AA1  0.5141629         0.7016303
#> 6:  Merged.M7     DNM2  0.5718974         0.6381931
```
