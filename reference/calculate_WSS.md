# Calculate Weighted Sum Scores Using Signed Protein Values

This function calculates module scores by sign-flipping proteins based
on their biomarker associations (mean_beta) and weighting by alpha
centrality. This approach aligns proteins in the direction of core
biomarker progression and weights by network importance and biomarker
association strength.

## Usage

``` r
calculate_WSS(
  expr_matrix,
  result = NULL,
  prebuilt = "core_AD_plasma_biomarkers"
)
```

## Arguments

- expr_matrix:

  A numeric matrix with samples as rows and genes/proteins as columns
  (INT-normalized). Column names should match the 'symbol' values in the
  result table.

- result:

  A data.table containing module information with required columns:
  module, symbol (gene identifiers), mean_beta (for sign-flipping),
  mean_alpha_scaled (for weighting). If `NULL` (the default), a prebuilt
  reference set shipped with the package is used instead; see
  `prebuilt`.

- prebuilt:

  Name of a prebuilt reference set to use when `result` is `NULL`.
  Defaults to `"core_AD_plasma_biomarkers"`. See
  [`list_prebuilt_wss()`](https://Longo-Lab.github.io/WSSmod/reference/list_prebuilt_wss.md)
  for available options. Ignored if `result` is supplied.

## Value

A list containing:

- scores:

  Matrix of weighted module scores (samples x modules)

- summary_stats:

  data.table with summary statistics for each module score

- module_sd:

  Named vector of standard deviations for each module

- module_proteins:

  List of protein names in each module

## Details

For each module, proteins are sign-flipped based on sign(mean_beta) to
align all proteins in the direction of disease progression. The signed
values are then weighted by mean_alpha_scaled (network centrality +
association strength) and averaged to produce a module score.

Higher module scores indicate greater module activity in the direction
of disease progression.

## See also

[`list_prebuilt_wss()`](https://Longo-Lab.github.io/WSSmod/reference/list_prebuilt_wss.md),
[`load_prebuilt_wss()`](https://Longo-Lab.github.io/WSSmod/reference/load_prebuilt_wss.md),
[`wss_prebuilt_terms()`](https://Longo-Lab.github.io/WSSmod/reference/wss_prebuilt_terms.md)

## Examples

``` r
expr_matrix <- matrix(
  c(1, 2, 3, 4, 5, 6, 7, 8),
  nrow = 2,
  dimnames = list(c("sample1", "sample2"), c("geneA", "geneB", "geneC", "geneD"))
)
cluster_result <- data.table::data.table(
  module = c("mod1", "mod1", "mod1", "mod2"),
  symbol = c("geneA", "geneB", "geneC", "geneD"),
  mean_beta = c(1, -1, 1, 1),
  mean_alpha_scaled = c(1, 1, 2, 1)
)

# Calculate weighted module scores
module_results <- calculate_WSS(expr_matrix, cluster_result)

# Access score matrix
scores <- module_results$scores

# View summary statistics
print(module_results$summary_stats)
#>    module   Min    Q1 Median  Mean    Q3   Max
#>    <char> <num> <num>  <num> <num> <num> <num>
#> 1:   mod1     2 2.125   2.25  2.25 2.375   2.5
#> 2:   mod2     7 7.250   7.50  7.50 7.750   8.0

# Check standard deviations
print(module_results$module_sd)
#>      mod1      mod2 
#> 0.3535534 0.7071068 
```
