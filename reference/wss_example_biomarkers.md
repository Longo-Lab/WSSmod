# Simulated Plasma Biomarkers for WSS Prediction Examples

A small, fully simulated table of raw plasma biomarker values for trying
out
[`normalize_wss_biomarkers()`](https://Longo-Lab.github.io/WSSmod/reference/normalize_wss_biomarkers.md)
and
[`predict_WSS()`](https://Longo-Lab.github.io/WSSmod/reference/predict_WSS.md).
Biomarkers are drawn from log-normal distributions at plausible scales.
The data are random draws and contain no participant data; they carry no
biological meaning and are for demonstration only.

## Usage

``` r
wss_example_biomarkers
```

## Format

A data.frame with 20 rows (named `patient01` to `patient20`) and 10
columns:

- Age:

  Age in years.

- Gender:

  `1` = male, `0` = not male.

- PlasmaPTau181, PlasmaAB142P, PlasmaAB140P, PlasmaABRatio,
  PlasmapTau217, PlasmapTau217_AB42Ratio, PlasmaGFAP, PlasmaNfL:

  Raw (un-normalized) plasma biomarker values.

## Source

Simulated; see `data-raw/simulate_examples.R` in the package source.

## See also

[`normalize_wss_biomarkers()`](https://Longo-Lab.github.io/WSSmod/reference/normalize_wss_biomarkers.md),
[`predict_WSS()`](https://Longo-Lab.github.io/WSSmod/reference/predict_WSS.md),
[wss_example_expression](https://Longo-Lab.github.io/WSSmod/reference/wss_example_expression.md)

## Examples

``` r
head(wss_example_biomarkers)
#>           Age Gender PlasmaPTau181 PlasmaAB142P PlasmaAB140P PlasmaABRatio
#> patient01  77      1        1.4203      28.6814     249.8937        0.0963
#> patient02  82      1        2.4113      25.5083     185.1647        0.1052
#> patient03  58      1        2.0431      24.1875     516.0148        0.0783
#> patient04  57      1        1.0156      21.9691     306.6212        0.1257
#> patient05  57      1        1.4921      24.8697     296.5344        0.0672
#> patient06  79      0        1.5639      31.7899     253.7137        0.0976
#>           PlasmapTau217 PlasmapTau217_AB42Ratio PlasmaGFAP PlasmaNfL
#> patient01        0.2373                  0.0053    53.2866   27.6111
#> patient02        0.1343                  0.0017    13.7432   34.9938
#> patient03        0.5072                  0.0175    38.1114   40.8202
#> patient04        0.1265                  0.0237   103.4999   26.0521
#> patient05        0.3000                  0.0276   123.2823   16.6924
#> patient06        0.0547                  0.0032    69.1005   13.4218
```
