## Test environments

* Local: Ubuntu 22.04 (WSL2), R 4.4.3
* GitHub Actions (`R-CMD-check`): macOS (release), Windows (release), Ubuntu (devel, release, oldrel-1)
* R-hub: Ubuntu, Windows and macOS (R-devel)
* win-builder: R-devel (2026-10-21 r90733) and R 4.6.1 (release)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.

The win-builder checks also report possibly misspelled words in DESCRIPTION:
"WSS" (the abbreviation for Weighted Sum Score), "proteomic" (a standard
term), and "et al" (in the method reference, Butler et al. (2026)
<doi:10.64898/2026.08.14.26360394>). All are intended.

(Two further messages appear only in the local sandbox and are unrelated to
the package: the note "unable to verify current time", and a warning that
`qpdf` is not installed.)

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
