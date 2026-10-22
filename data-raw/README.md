# data-raw

Scripts recording how the shipped data were made. This directory is excluded
from the built package (`.Rbuildignore`).

| Script | Produces | Runnable by |
|--------|----------|-------------|
| `simulate_examples.R` | `data/wss_example_expression.rda`, `data/wss_example_biomarkers.rda` (fully simulated, seed 20261001) | anyone, from the package root |
| `extdata_core_AD_plasma_biomarkers.R` | the four `inst/extdata/core_AD_plasma_biomarkers.*` files | Longo Lab only (needs the private upstream Box data) |

`extdata_core_AD_plasma_biomarkers.R` reproduces all four `inst/extdata` files
byte-for-byte (checked by MD5 against the shipped files). The two module tables
are unchanged copies from the paper pipeline, the model is the
`Biomarkers_Age_Sex` variant extracted from the upstream final models, and the
biomarker reference is the 389-row joint complete-case baseline set with each
biomarker independently shuffled (seed 42). The script is documentation of
provenance for the private data; the participant-level inputs are not
distributed.
