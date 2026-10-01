## Test environments

* Local: Ubuntu 22.04 (WSL2), R 4.4.3
* GitHub Actions (`R-CMD-check`): macOS (release), Windows (release), Ubuntu (devel, release, oldrel-1)

## R CMD check results

`R CMD check --as-cran` (local, R 4.4.3):

0 errors | 0 warnings | 1 note

* This is a new submission.

(A second note, "unable to verify current time", appears only in the local
sandbox and is unrelated to the package.)

## Notes for reviewers

* `glmnet`, `joinet` and `RNOmni` are in Suggests. They are only needed for
  `predict_WSS()` and `normalize_wss_biomarkers(method = "self")`. Code,
  examples and tests that use them are guarded with `requireNamespace()` /
  `skip_if_not_installed()`, and the model-dependent tests also use
  `skip_on_cran()`.
* The package ships small data files in `inst/extdata` (a module weighting
  table, module term labels, a fitted prediction model, and an anonymized
  baseline biomarker reference) and two simulated example datasets in `data/`.
  The example datasets contain no participant data.

## Downstream dependencies

There are no reverse dependencies (new package).
