# Simulated Expression Matrix for WSS Examples

A small, fully simulated samples-by-proteins matrix for trying out
[`calculate_WSS()`](https://Longo-Lab.github.io/WSSmod/reference/calculate_WSS.md).
Its columns are the 887 proteins (gene symbols) that survive the default
deduplication and module-size filter of the bundled
`"core_AD_plasma_biomarkers"` reference set, so
`calculate_WSS(wss_example_expression)` scores all 75 modules without
warnings. Values are on a z-score-like scale, as expected for normalized
input.

## Usage

``` r
wss_example_expression
```

## Format

A numeric matrix with 50 rows (samples, named `sample01` to `sample50`)
and 887 columns (proteins, named by gene symbol).

## Source

Simulated; see `data-raw/simulate_examples.R` in the package source.

## Details

Each module was given a latent per-sample factor, and its proteins load
on that factor in the direction of their `mean_beta`, plus independent
noise. The data are random draws and contain no participant data; they
carry no biological meaning and are for demonstration only.

## See also

[`calculate_WSS()`](https://Longo-Lab.github.io/WSSmod/reference/calculate_WSS.md),
[wss_example_biomarkers](https://Longo-Lab.github.io/WSSmod/reference/wss_example_biomarkers.md)

## Examples

``` r
dim(wss_example_expression)
#> [1]  50 887
res <- calculate_WSS(wss_example_expression)
head(res$scores[, 1:3])
#>          Merged.M17   Merged.M7   Merged.M4
#> sample01 -0.8723738  0.36490845  0.04891509
#> sample02  0.2344854 -0.50295921 -0.54133339
#> sample03  0.8852222 -0.41936437 -0.56949082
#> sample04  0.8491436  0.05028353 -0.97018547
#> sample05 -0.2620853 -1.56467707  0.53368294
#> sample06 -0.6801977 -1.26890805  0.48510968
```
