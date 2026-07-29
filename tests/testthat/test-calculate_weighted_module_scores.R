test_that("calculate_weighted_module_scores computes expected weighted averages", {
  set.seed(1)
  expr_matrix <- matrix(
    c(1, 2, 3, 4, 5, 6, 7, 8),
    nrow = 2,
    dimnames = list(c("sample1", "sample2"), c("geneA", "geneB", "geneC", "geneD"))
  )

  result <- data.table::data.table(
    module = c("mod1", "mod1", "mod1", "mod2"),
    symbol = c("geneA", "geneB", "geneC", "geneD"),
    mean_beta = c(1, -1, 1, 1),
    mean_alpha_scaled = c(1, 1, 2, 1)
  )

  out <- calculate_weighted_module_scores(expr_matrix, result)

  expect_named(out, c("scores", "summary_stats", "module_sd", "module_proteins"))
  expect_equal(colnames(out$scores), c("mod1", "mod2"))

  # sample1: geneA=1, geneB=3 -> -3 (flipped), geneC=5, weights 1,1,2 -> (1 - 3 + 5*2)/4 = 2
  expect_equal(out$scores["sample1", "mod1"], 2)
  expect_equal(out$scores["sample1", "mod2"], 7)

  expect_equal(out$module_proteins$mod1, c("geneA", "geneB", "geneC"))
  expect_equal(out$module_proteins$mod2, "geneD")
})

test_that("calculate_weighted_module_scores falls back to equal weights when all weights are zero", {
  expr_matrix <- matrix(
    c(1, 2, 3, 4),
    nrow = 2,
    dimnames = list(c("sample1", "sample2"), c("geneA", "geneB"))
  )

  result <- data.table::data.table(
    module = c("mod1", "mod1"),
    symbol = c("geneA", "geneB"),
    mean_beta = c(1, 1),
    mean_alpha_scaled = c(0, 0)
  )

  expect_warning(
    out <- calculate_weighted_module_scores(expr_matrix, result),
    "missing or zero weights"
  )

  expect_equal(out$scores["sample1", "mod1"], mean(c(1, 3)))
})

test_that("calculate_weighted_module_scores warns when a module has no matching proteins", {
  expr_matrix <- matrix(
    c(1, 2),
    nrow = 2,
    dimnames = list(c("sample1", "sample2"), "geneA")
  )

  result <- data.table::data.table(
    module = "mod1",
    symbol = "geneMissing",
    mean_beta = 1,
    mean_alpha_scaled = 1
  )

  expect_warning(
    calculate_weighted_module_scores(expr_matrix, result),
    "no proteins in expression matrix"
  )
})
