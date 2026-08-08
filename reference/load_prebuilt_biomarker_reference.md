# Load a Prebuilt Biomarker Reference Distribution

Some prebuilt reference sets ship a baseline biomarker reference
distribution used by
[`normalize_wss_biomarkers()`](https://Longo-Lab.github.io/WSSmod/reference/normalize_wss_biomarkers.md)
(`method = "project"`) to rank-normalize new raw biomarker values
without needing a sizeable cohort of your own. Each biomarker's
reference values are drawn from the real baseline cohort, but
**independently randomly shuffled** so that no two values across
different biomarkers can be attributed to the same original subject.
This preserves each biomarker's marginal distribution exactly (every
value, tie, and extreme is real) while destroying the joint/multivariate
fingerprint that carries re-identification risk.

## Usage

``` r
load_prebuilt_biomarker_reference(prebuilt = "core_AD_plasma_biomarkers")
```

## Arguments

- prebuilt:

  Name of a prebuilt reference set. See
  [`list_prebuilt_wss()`](https://Longo-Lab.github.io/WSSmod/reference/list_prebuilt_wss.md)
  for available options.

## Value

A data.table with columns `biomarker` and `value` (long format, one row
per reference observation; biomarkers have differing numbers of rows due
to differing missingness in the original cohort), or `NULL` if the given
reference set has no associated biomarker reference.

## Details

Because of this, **the returned reference has no valid joint structure**
— inter-biomarker correlations are deliberately destroyed. Only use it
for independent per-biomarker operations (rank/quantile lookups), never
for anything requiring multivariate relationships (e.g. correlation or
joint modeling).

## See also

[`normalize_wss_biomarkers()`](https://Longo-Lab.github.io/WSSmod/reference/normalize_wss_biomarkers.md),
[`project_rank_norm()`](https://Longo-Lab.github.io/WSSmod/reference/project_rank_norm.md)
