# Get Module Term Labels for a Prebuilt WSS Reference Set

Some prebuilt reference sets ship a companion table of human-readable
term labels for each module, summarized from the underlying pathway
enrichment used to build the module.

## Usage

``` r
wss_prebuilt_terms(prebuilt = "core_AD_plasma_biomarkers")
```

## Arguments

- prebuilt:

  Name of a prebuilt reference set. See
  [`list_prebuilt_wss()`](https://Longo-Lab.github.io/WSSmod/reference/list_prebuilt_wss.md)
  for available options.

## Value

A data.table of module term labels, or `NULL` if the given reference set
has no associated terms file.

## Examples

``` r
head(wss_prebuilt_terms())
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
