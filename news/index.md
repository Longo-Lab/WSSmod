# Changelog

## WSSmod 0.1.0

Initial versioned release.

### Module scores

- [`calculate_WSS()`](https://Longo-Lab.github.io/WSSmod/reference/calculate_WSS.md)
  computes Weighted Sum of Scores module-level summaries from proteomic
  data, using STRING network connectivity and biomarker-association
  weighting to align and aggregate protein-level signal into module
  scores. Supports `dedup_by` and `min_module_size` options for
  controlling how genes assigned to multiple modules are deduplicated
  and how small modules are filtered.
- [`load_prebuilt_wss()`](https://Longo-Lab.github.io/WSSmod/reference/load_prebuilt_wss.md),
  [`list_prebuilt_wss()`](https://Longo-Lab.github.io/WSSmod/reference/list_prebuilt_wss.md),
  and
  [`wss_prebuilt_terms()`](https://Longo-Lab.github.io/WSSmod/reference/wss_prebuilt_terms.md)
  provide access to a bundled reference module/weighting table
  (`"core_AD_plasma_biomarkers"`) and its module term labels.

### Prediction from core biomarkers

- [`predict_WSS()`](https://Longo-Lab.github.io/WSSmod/reference/predict_WSS.md)
  predicts all 75 module scores for the `"core_AD_plasma_biomarkers"`
  reference set from Age, Gender, and 8 core plasma biomarkers, using a
  bundled `joinet` two-layer elastic net model, without requiring the
  full proteomic panel.
- [`normalize_wss_biomarkers()`](https://Longo-Lab.github.io/WSSmod/reference/normalize_wss_biomarkers.md)
  builds the rank-normalized biomarker inputs
  [`predict_WSS()`](https://Longo-Lab.github.io/WSSmod/reference/predict_WSS.md)
  requires from raw biomarker values, via either `method = "project"`
  (projecting onto a bundled, anonymized baseline reference
  distribution, using
  [`project_rank_norm()`](https://Longo-Lab.github.io/WSSmod/reference/project_rank_norm.md)
  – supports a single new patient) or `method = "self"`
  (rank-normalizing within your own cohort).
- [`project_rank_norm()`](https://Longo-Lab.github.io/WSSmod/reference/project_rank_norm.md)
  and
  [`load_prebuilt_biomarker_reference()`](https://Longo-Lab.github.io/WSSmod/reference/load_prebuilt_biomarker_reference.md)
  expose the underlying projection utility and reference distribution
  directly.

The bundled model and reference data were verified end-to-end against
the original analysis’s held-out follow-up cohort, reproducing published
predictions to floating-point precision.
