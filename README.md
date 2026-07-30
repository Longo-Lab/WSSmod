
<!-- README.md is generated from README.Rmd. Please edit that file -->

# WSSmod

<!-- badges: start -->

[![Version](https://img.shields.io/badge/version-0.0.0.9000-blue.svg)](https://github.com/Longo-Lab/WSSmod/blob/main/DESCRIPTION)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/Longo-Lab/WSSmod/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Longo-Lab/WSSmod/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/Longo-Lab/WSSmod/branch/main/graph/badge.svg)](https://app.codecov.io/gh/Longo-Lab/WSSmod?branch=main)
<!-- badges: end -->

WSSmod provides tools to construct Weighted Sum of Scores (WSS)
module-level summaries from proteomic data. Proteins are grouped into
STRING network modules, sign-flipped to align with the direction of core
biomarker association, and combined into a single module score per
sample using a weighted average based on network connectivity (alpha
centrality) and biomarker association strength.

This package accompanies the Longo Lab plasma biomarker proteomics work.

## Installation

You can install the development version of WSSmod from GitHub with:

``` r
# install.packages("pak")
pak::pak("Longo-Lab/WSSmod")
```

## Example

Using a custom module/weighting table:

``` r
library(WSSmod)

expr_matrix <- matrix(
  c(1, 2, 3, 4, 5, 6, 7, 8),
  nrow = 2,
  dimnames = list(c("sample1", "sample2"), c("geneA", "geneB", "geneC", "geneD"))
)

result <- data.table::data.table(
  module = c("mod1", "mod1", "mod1", "mod2"),
  symbol = c("geneA", "geneB", "geneC", "geneD"),
  mean_beta = c(1, -1, 1, 1),
  mean_alpha_scaled = c(1, 1, 2, 1)
)

WSSmod::calculate_WSS(expr_matrix, result, min_module_size = 1)
#> $scores
#>         mod1 mod2
#> sample1  2.0    7
#> sample2  2.5    8
#> 
#> $summary_stats
#>    module   Min    Q1 Median  Mean    Q3   Max
#>    <char> <num> <num>  <num> <num> <num> <num>
#> 1:   mod1     2 2.125   2.25  2.25 2.375   2.5
#> 2:   mod2     7 7.250   7.50  7.50 7.750   8.0
#> 
#> $module_sd
#>      mod1      mod2 
#> 0.3535534 0.7071068 
#> 
#> $module_proteins
#> $module_proteins$mod1
#> [1] "geneC" "geneA" "geneB"
#> 
#> $module_proteins$mod2
#> [1] "geneD"
```

WSSmod also has prebuilt module/weighting reference sets, including those from the Butler et al paper, so you don’t
have to build your own module table from scratch. `calculate_WSS()` uses
`prebuilt = "core_AD_plasma_biomarkers"` by default when `result =
NULL`:

``` r
list_prebuilt_wss()

# expr_matrix here has samples as rows and protein/gene symbols as columns,
# matching the symbols used by the prebuilt reference set
WSSmod::calculate_WSS(expr_matrix)

# Human-readable module term labels for the same reference set
wss_prebuilt_terms("core_AD_plasma_biomarkers")
```
