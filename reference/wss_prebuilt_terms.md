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
