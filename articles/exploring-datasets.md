# Exploring the bundled datasets

``` r

library(WSSmod)
```

WSSmod ships files in `inst/extdata` that belong to a named *prebuilt
reference set*, plus two small simulated datasets for examples. This
vignette describes each one: what it is, how to load it, and what it can
and cannot be used for. The other vignettes put them to use:
[Calculating Weighted Sum
Scores](https://Longo-Lab.github.io/WSSmod/articles/calculating-wss.md)
and [Predicting Weighted Sum Scores from plasma
biomarkers](https://Longo-Lab.github.io/WSSmod/articles/predicting-wss.md).

## Available reference sets

``` r

list_prebuilt_wss()
#> [1] "core_AD_plasma_biomarkers"
```

There is one set, `"core_AD_plasma_biomarkers"`. Every function below
takes it as its `prebuilt` argument, and as the default. A set can have
up to four files, and each is read with its own function:

| File | Function | Contents |
|----|----|----|
| scores | [`load_prebuilt_wss()`](https://Longo-Lab.github.io/WSSmod/reference/load_prebuilt_wss.md) | module and protein weighting table |
| terms | [`wss_prebuilt_terms()`](https://Longo-Lab.github.io/WSSmod/reference/wss_prebuilt_terms.md) | module labels and enrichment summary |
| model | [`load_prebuilt_wss_model()`](https://Longo-Lab.github.io/WSSmod/reference/load_prebuilt_wss_model.md) | fitted `joinet` prediction model |
| biomarker reference | [`load_prebuilt_biomarker_reference()`](https://Longo-Lab.github.io/WSSmod/reference/load_prebuilt_biomarker_reference.md) | anonymized baseline biomarker values |

## The module and protein weighting table

[`load_prebuilt_wss()`](https://Longo-Lab.github.io/WSSmod/reference/load_prebuilt_wss.md)
returns the table that
[`calculate_WSS()`](https://Longo-Lab.github.io/WSSmod/reference/calculate_WSS.md)
uses by default, with one row per protein-in-module assignment.

``` r

weights <- load_prebuilt_wss()
dim(weights)
#> [1] 992  13
names(weights)
#>  [1] "Biomarker"         "symbol"            "cluster"          
#>  [4] "GENE"              "ix"                "term"             
#>  [7] "orig_groups"       "orig_terms"        "mean_beta"        
#> [10] "mean_degree"       "mean_strength"     "mean_alpha_scaled"
#> [13] "module"
```

The columns
[`calculate_WSS()`](https://Longo-Lab.github.io/WSSmod/reference/calculate_WSS.md)
needs are

- `module`: module name, built as `Biomarker`, a dot, and `cluster` (for
  example `Merged.M17`),
- `symbol`: protein gene symbol, matching the columns of your matrix,
- `mean_beta`: its sign sets the direction a protein is aligned in,
- `mean_alpha_scaled`: the weight of the protein within its module,
  scaled to at most 1.

``` r

head(weights[, c("module", "symbol", "term", "mean_beta", "mean_alpha_scaled")])
#>        module   symbol                      term  mean_beta mean_alpha_scaled
#>        <char>   <char>                    <char>      <num>             <num>
#> 1: Merged.M17   ZWILCH                Cell cycle -0.4732375         1.0000000
#> 2:  Merged.M7    HYOU1 Chaperone|Protein folding  0.1603553         0.9961153
#> 3:  Merged.M4    EFNB2          Ephrin signaling  0.4591349         0.8693081
#> 4:  Merged.M7  DNAJB11 Chaperone|Protein folding  0.2211760         0.8307261
#> 5:  Merged.M7 HSP90AA1 Chaperone|Protein folding  0.5141629         0.7016303
#> 6:  Merged.M7     DNM2 Chaperone|Protein folding  0.5718974         0.6381931
```

The other columns describe where a row came from: the identifiers `GENE`
and `ix` (an aptamer identifier), the module’s `term` label, the
original biomarker-specific groups and terms that the protein came from
(`orig_groups` and `orig_terms`), and network connectivity summaries
(`mean_degree` and `mean_strength`).

### Why a protein can be in several modules

A protein can be assigned to more than one module, because the modules
were built from several biomarker analyses and then merged.

``` r

c(
  rows = nrow(weights),
  proteins = length(unique(weights$symbol)),
  modules = length(unique(weights$module))
)
#>     rows proteins  modules 
#>      992      898       82

# Proteins that are in more than one module
multi <- weights[, .N, by = symbol][N > 1]
nrow(multi)
#> [1] 89
weights[symbol == multi$symbol[1], c("module", "symbol", "term", "mean_alpha_scaled")]
#>       module symbol                             term mean_alpha_scaled
#>       <char> <char>                           <char>             <num>
#> 1: Merged.M1  EPHA2 PI3K-Akt|JAK-STAT|MAPK signaling         0.5786472
#> 2: Merged.M4  EPHA2                 Ephrin signaling         0.2760032
```

[`calculate_WSS()`](https://Longo-Lab.github.io/WSSmod/reference/calculate_WSS.md)
removes this overlap before scoring. It keeps each protein in the module
where it has the highest `mean_alpha_scaled` (the `dedup_by` argument),
then drops modules with fewer than `min_module_size` proteins (3 by
default). That leaves the 75 modules that the prediction model also
covers.

``` r

res <- calculate_WSS(wss_example_expression)
c(
  modules_scored = ncol(res$scores),
  proteins_scored = length(unique(unlist(res$module_proteins)))
)
#>  modules_scored proteins_scored 
#>              75             887
```

## Module term labels

[`wss_prebuilt_terms()`](https://Longo-Lab.github.io/WSSmod/reference/wss_prebuilt_terms.md)
returns one row per scored module with a short label, the group the
module came from, and `mean.Zsum`, the Zsummary preservation score of
the module across 12 independent GNPC cohorts in the
[preprint](https://www.medrxiv.org/content/10.64898/2026.08.14.26360394v1.full-text).

``` r

terms <- wss_prebuilt_terms()
dim(terms)
#> [1] 75  4
head(terms)
#>                name                                        term         group
#>              <char>                                      <char>        <char>
#> 1:        Merged.M1            PI3K-Akt|JAK-STAT|MAPK signaling        Merged
#> 2:       Merged.M11                          ER-Golgi transport        Merged
#> 3:        Merged.M7                   Chaperone|Protein folding        Merged
#> 4:        Merged.M9                           Carbon metabolism        Merged
#> 5:        Merged.M2 Adaptive immunity|Immunoglobulins|Cytokines        Merged
#> 6: PlasmaPTau181.M8                      Protein Ubiquitination PlasmaPTau181
#>    mean.Zsum
#>        <num>
#> 1: 10.392179
#> 2:  8.924015
#> 3:  8.463061
#> 4:  7.078961
#> 5:  6.894285
#> 6:  6.489790
```

The `name` column matches the module names in the weighting table and in
the columns of
[`calculate_WSS()`](https://Longo-Lab.github.io/WSSmod/reference/calculate_WSS.md)
and
[`predict_WSS()`](https://Longo-Lab.github.io/WSSmod/reference/predict_WSS.md)
results, so it can be used to label them:

``` r

modules <- colnames(res$scores)[1:3]
terms$term[match(modules, terms$name)]
#> [1] "Cell cycle"                "Chaperone|Protein folding"
#> [3] "Ephrin signaling"
```

`group` is the biomarker analysis a module came from, with `Merged` for
modules built from all of them together.

``` r

table(terms$group)
#> 
#>        Merged  PlasmaAB140P  PlasmaAB142P    PlasmaGFAP     PlasmaNfL 
#>            32            15             4             2             5 
#> PlasmaPTau181 PlasmapTau217 
#>            16             1
```

## The prediction model

[`load_prebuilt_wss_model()`](https://Longo-Lab.github.io/WSSmod/reference/load_prebuilt_wss_model.md)
returns a list with the fitted model and what it needs.

``` r

model_info <- load_prebuilt_wss_model()
names(model_info)
#> [1] "model"      "outcomes"   "variant"    "alpha_base" "alpha_meta"
#> [6] "x_cols"
model_info$variant
#> [1] "Biomarkers_Age_Sex"
model_info$x_cols
#>  [1] "Age"                     "Gender"                 
#>  [3] "PlasmaPTau181"           "PlasmaAB142P"           
#>  [5] "PlasmaAB140P"            "PlasmaABRatio"          
#>  [7] "PlasmapTau217"           "PlasmapTau217_AB42Ratio"
#>  [9] "PlasmaGFAP"              "PlasmaNfL"
head(model_info$outcomes)
#> [1] "Merged.M17"        "Merged.M7"         "Merged.M4"        
#> [4] "PlasmaPTau181.M26" "Merged.M8"         "Merged.M1"
length(model_info$outcomes)
#> [1] 75
```

- `x_cols` are the predictor columns
  [`predict_WSS()`](https://Longo-Lab.github.io/WSSmod/reference/predict_WSS.md)
  requires, in order: age, sex and the eight plasma biomarkers.
- `outcomes` are the 75 module names it predicts, the same names as in
  the terms table.
- `variant` names the set of predictors the model was trained on.
- `model` is the fitted `joinet` object. Predicting from it needs the
  suggested packages `glmnet` and `joinet`, which is why the other
  vignette guards those chunks. Reading the list does not need them, so
  nothing above does.

To predict, use
[`predict_WSS()`](https://Longo-Lab.github.io/WSSmod/reference/predict_WSS.md)
rather than the model directly.

## The anonymized biomarker reference

[`load_prebuilt_biomarker_reference()`](https://Longo-Lab.github.io/WSSmod/reference/load_prebuilt_biomarker_reference.md)
returns baseline values for the eight plasma biomarkers, in long format,
that
[`normalize_wss_biomarkers()`](https://Longo-Lab.github.io/WSSmod/reference/normalize_wss_biomarkers.md)
uses with `method = "project"`.

``` r

reference <- load_prebuilt_biomarker_reference()
head(reference, 3)
#>        biomarker value
#>           <char> <num>
#> 1: PlasmaPTau181  0.89
#> 2: PlasmaPTau181  1.35
#> 3: PlasmaPTau181  1.48
reference[, .(n = .N, median = signif(median(value), 3)), by = biomarker]
#>                  biomarker     n   median
#>                     <char> <int>    <num>
#> 1:           PlasmaPTau181   389 1.62e+00
#> 2:            PlasmaAB142P   389 2.78e+01
#> 3:            PlasmaAB140P   389 3.03e+02
#> 4:           PlasmaABRatio   389 9.22e-02
#> 5:           PlasmapTau217   389 1.49e-01
#> 6: PlasmapTau217_AB42Ratio   389 5.14e-03
#> 7:              PlasmaGFAP   389 5.88e+01
#> 8:               PlasmaNfL   389 2.03e+01
```

The values are real, but each biomarker’s values were shuffled
independently, so the reference keeps every biomarker’s distribution
exactly while removing the link between biomarkers that could identify a
person. Use it only for per-biomarker questions, such as where a new
value falls within the distribution. It has no valid joint structure, so
do not use it to look at relationships between biomarkers or to train or
evaluate models.

## Simulated example data

Two datasets, available as soon as the package is loaded, are simulated
and have no biological meaning. They exist for examples and tests:

- `wss_example_expression`: a 50-sample by 887-protein matrix covering
  every protein in the default modules,
- `wss_example_biomarkers`: age, sex and eight raw biomarkers for 20
  patients.

``` r

dim(wss_example_expression)
#> [1]  50 887
dim(wss_example_biomarkers)
#> [1] 20 10
```

## Provenance and use

The prebuilt files were derived from plasma proteomic and biomarker data
from the Longo Lab, as described in the
[preprint](https://www.medrxiv.org/content/10.64898/2026.08.14.26360394v1.full-text).
They are included so that scores can be calculated and predicted without
access to the original data. (Maintainer to finalize the provenance and
licensing statement for each file.)

## Adding a reference set

Reference sets are listed in one registry in the package source,
`.wss_prebuilt_registry` in `R/prebuilt_wss.R`. A new set is registered
by giving it a name and the file names of its parts. Only the scores
file is required; the terms, model and biomarker reference files are
optional, and the loaders return `NULL` for a part a set does not have.
The files go in `inst/extdata`, and the new name then appears in
[`list_prebuilt_wss()`](https://Longo-Lab.github.io/WSSmod/reference/list_prebuilt_wss.md).
