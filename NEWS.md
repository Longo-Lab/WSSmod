# WSSmod 0.1.0

Initial versioned release.

## Module scores

* `calculate_WSS()` computes Weighted Sum of Scores module-level summaries
  from proteomic data, using STRING network connectivity and
  biomarker-association weighting to align and aggregate protein-level
  signal into module scores. Supports `dedup_by` and `min_module_size`
  options for controlling how genes assigned to multiple modules are
  deduplicated and how small modules are filtered.
* `load_prebuilt_wss()`, `list_prebuilt_wss()`, and `wss_prebuilt_terms()`
  provide access to a bundled reference module/weighting table
  (`"core_AD_plasma_biomarkers"`) and its module term labels.

## Prediction from core biomarkers

* `predict_WSS()` predicts all 75 module scores for the
  `"core_AD_plasma_biomarkers"` reference set from Age, Gender, and 8 core
  plasma biomarkers, using a bundled `joinet` two-layer elastic net model,
  without requiring the full proteomic panel.
* `normalize_wss_biomarkers()` builds the rank-normalized biomarker inputs
  `predict_WSS()` requires from raw biomarker values, via either
  `method = "project"` (projecting onto a bundled, anonymized baseline
  reference distribution, using `project_rank_norm()` -- supports a single
  new patient) or `method = "self"` (rank-normalizing within your own
  cohort).
* `project_rank_norm()` and `load_prebuilt_biomarker_reference()` expose
  the underlying projection utility and reference distribution directly.

The bundled model and reference data were verified end-to-end against the
original analysis's held-out follow-up cohort, reproducing published
predictions to floating-point precision.
