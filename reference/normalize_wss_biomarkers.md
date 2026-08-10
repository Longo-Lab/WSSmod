# Normalize Raw Biomarker Values for predict_WSS()

Builds the rank-based inverse-normal-transformed biomarker columns
[`predict_WSS()`](https://Longo-Lab.github.io/WSSmod/reference/predict_WSS.md)
requires, from raw biomarker values, using one of two methods. Output
columns keep the same names as the input columns (e.g. `PlasmaPTau181`
stays `PlasmaPTau181`) –
[`predict_WSS()`](https://Longo-Lab.github.io/WSSmod/reference/predict_WSS.md)
reads whatever column names its underlying model actually expects.

## Usage

``` r
normalize_wss_biomarkers(
  raw_biomarkers,
  prebuilt = "core_AD_plasma_biomarkers",
  method = c("project", "self")
)
```

## Arguments

- raw_biomarkers:

  A data.frame or matrix, one row per sample, with columns named by the
  raw (unsuffixed) biomarker names (e.g. `PlasmaPTau181`, not
  `PlasmaPTau181_norm`).

- prebuilt:

  Name of a prebuilt reference set to use when `method = "project"`. See
  [`list_prebuilt_wss()`](https://Longo-Lab.github.io/WSSmod/reference/list_prebuilt_wss.md)
  for available options. Ignored when `method = "self"`.

- method:

  Normalization method:

  `"project"`

  :   (default) Project each raw value onto a bundled reference
      distribution via
      [`project_rank_norm()`](https://Longo-Lab.github.io/WSSmod/reference/project_rank_norm.md).
      Works for a single sample or a small/differently-distributed
      cohort, since each value is scored independently against the fixed
      reference.

  `"self"`

  :   Rank-normalize (via
      [`RNOmni::RankNorm()`](https://rdrr.io/pkg/RNOmni/man/RankNorm.html))
      within `raw_biomarkers` itself. Appropriate if you have a sizeable
      cohort of your own with a distribution you're comfortable
      normalizing against directly, instead of the bundled reference.
      Requires at least two samples, since a rank transform needs
      multiple values to rank against.

## Value

A data.frame with one normalized column per input biomarker column (same
names as `raw_biomarkers`), row names preserved. Combine with
`Age`/`Gender` columns and pass to
[`predict_WSS()`](https://Longo-Lab.github.io/WSSmod/reference/predict_WSS.md).

## See also

[`predict_WSS()`](https://Longo-Lab.github.io/WSSmod/reference/predict_WSS.md),
[`project_rank_norm()`](https://Longo-Lab.github.io/WSSmod/reference/project_rank_norm.md),
[`load_prebuilt_biomarker_reference()`](https://Longo-Lab.github.io/WSSmod/reference/load_prebuilt_biomarker_reference.md)

## Examples

``` r
raw_biomarkers <- data.frame(PlasmaPTau181 = 1.5, PlasmaNfL = 20)
normalize_wss_biomarkers(raw_biomarkers)
#>   PlasmaPTau181    PlasmaNfL
#> 1    -0.2168629 -0.009634851
```
