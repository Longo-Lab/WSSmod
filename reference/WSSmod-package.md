# WSSmod: Weighted Sum Scores for Protein Network Modules

Computes Weighted Sum Score (WSS) module-level summaries from proteomic
data. Protein-level values are sign-aligned to their association with
core biomarkers and combined within 'STRING' network modules, weighted
by network connectivity and association strength. Bundled reference
modules and a fitted 'joinet' model also allow all module scores to be
predicted from age, sex and a small panel of plasma biomarkers, without
measuring the full proteomic panel.

## Details

WSSmod computes Weighted Sum Score (WSS) values for protein network
modules, and predicts them from a small panel of plasma biomarkers.

**Scoring modules from proteomic data**

- [`calculate_WSS()`](https://Longo-Lab.github.io/WSSmod/reference/calculate_WSS.md)
  scores each module from a normalized and scaled samples-by-proteins
  matrix, using the bundled reference modules or your own.

- [`project_rank_norm()`](https://Longo-Lab.github.io/WSSmod/reference/project_rank_norm.md)
  projects new values onto a frozen rank-normal reference, for
  normalizing new samples against a reference cohort.

**Predicting module scores from plasma biomarkers**

- [`normalize_wss_biomarkers()`](https://Longo-Lab.github.io/WSSmod/reference/normalize_wss_biomarkers.md)
  rank-normalizes raw biomarker values, for a single patient or a
  cohort.

- [`predict_WSS()`](https://Longo-Lab.github.io/WSSmod/reference/predict_WSS.md)
  predicts all module scores from age, sex and the normalized biomarkers
  (needs the suggested packages `glmnet` and `joinet`).

**Bundled reference data**

- [`list_prebuilt_wss()`](https://Longo-Lab.github.io/WSSmod/reference/list_prebuilt_wss.md),
  [`load_prebuilt_wss()`](https://Longo-Lab.github.io/WSSmod/reference/load_prebuilt_wss.md),
  [`wss_prebuilt_terms()`](https://Longo-Lab.github.io/WSSmod/reference/wss_prebuilt_terms.md),
  [`load_prebuilt_wss_model()`](https://Longo-Lab.github.io/WSSmod/reference/load_prebuilt_wss_model.md)
  and
  [`load_prebuilt_biomarker_reference()`](https://Longo-Lab.github.io/WSSmod/reference/load_prebuilt_biomarker_reference.md)
  give access to the prebuilt reference set.

- [wss_example_expression](https://Longo-Lab.github.io/WSSmod/reference/wss_example_expression.md)
  and
  [wss_example_biomarkers](https://Longo-Lab.github.io/WSSmod/reference/wss_example_biomarkers.md)
  are small simulated datasets for trying out the functions.

## See also

Useful links:

- <https://github.com/Longo-Lab/WSSmod>

- <https://Longo-Lab.github.io/WSSmod/>

- Report bugs at <https://github.com/Longo-Lab/WSSmod/issues>

## Author

**Maintainer**: Robert Butler <rrbutler@stanford.edu>
([ORCID](https://orcid.org/0000-0003-0200-8682))

Other contributors:

- WSSmod authors \[copyright holder\]
