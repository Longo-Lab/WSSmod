test_that("list_prebuilt_wss reports the built-in reference set", {
  expect_true("core_AD_plasma_biomarkers" %in% list_prebuilt_wss())
})

test_that("load_prebuilt_wss returns a filtered, deduplicated module table", {
  result <- load_prebuilt_wss("core_AD_plasma_biomarkers")

  expect_true(data.table::is.data.table(result))
  expect_true(all(c("module", "symbol", "mean_beta", "mean_alpha_scaled") %in% names(result)))

  # No protein should appear in more than one module (deduplicated by symbol)
  expect_equal(anyDuplicated(result$symbol), 0)

  # Singleton/doubleton modules should have been filtered out
  sizes <- result[, .N, by = module]
  expect_true(all(sizes$N > 2))
})

test_that("load_prebuilt_wss errors on an unknown prebuilt name", {
  expect_error(
    load_prebuilt_wss("not_a_real_set"),
    "Unknown prebuilt WSS reference set"
  )
})

test_that("wss_prebuilt_terms returns module term labels for the default set", {
  terms <- wss_prebuilt_terms("core_AD_plasma_biomarkers")

  expect_true(data.table::is.data.table(terms))
  expect_true(all(c("name", "term", "group", "mean.Zsum") %in% names(terms)))
})
