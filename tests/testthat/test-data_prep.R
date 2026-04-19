# tests/testthat/test-data_prep.R
# Tests for data preparation functions
#
# Shared fixture using a subset of example_se (real GTEx data)

library(SummarizedExperiment)
data(example_se)
se <- example_se[1:50, 1:30]

# ── top_variable_features() ───────────────────────────────────────────────────

test_that("top_variable_features errors on non-SummarizedExperiment", {
  expect_error(
    top_variable_features("not_a_se"),
    "se must be a SummarizedExperiment object"
  )
})

test_that("top_variable_features errors on missing assay", {
  expect_error(
    top_variable_features(se, assay_name = "tpm"),
    "assay 'tpm' not found in se"
  )
})

test_that("top_variable_features returns a SummarizedExperiment", {
  result <- top_variable_features(se, n = 20)
  expect_s4_class(result, "SummarizedExperiment")
})

test_that("top_variable_features returns exactly n features", {
  result <- top_variable_features(se, n = 20)
  expect_equal(nrow(result), 20)
})

test_that("top_variable_features preserves sample count", {
  result <- top_variable_features(se, n = 20)
  expect_equal(ncol(result), ncol(se))
})

test_that("top_variable_features returns most variable genes", {
  result_10  <- top_variable_features(se, n = 10)
  result_20  <- top_variable_features(se, n = 20)
  expect_true(all(rownames(result_10) %in% rownames(result_20)))
})

# ── select_random_samples() ───────────────────────────────────────────────────

test_that("select_random_samples errors on non-SummarizedExperiment", {
  expect_error(
    select_random_samples("not_a_se"),
    "se must be a SummarizedExperiment object"
  )
})

test_that("select_random_samples errors on invalid n", {
  expect_error(
    select_random_samples(se, n = -1),
    "n must be a single positive integer"
  )
})

test_that("select_random_samples returns a SummarizedExperiment", {
  set.seed(42)
  result <- select_random_samples(se, n = 10)
  expect_s4_class(result, "SummarizedExperiment")
})

test_that("select_random_samples returns exactly n samples", {
  set.seed(42)
  result <- select_random_samples(se, n = 10)
  expect_equal(ncol(result), 10)
})

test_that("select_random_samples preserves feature count", {
  set.seed(42)
  result <- select_random_samples(se, n = 10)
  expect_equal(nrow(result), nrow(se))
})

test_that("select_random_samples caps at available samples", {
  set.seed(42)
  result <- select_random_samples(se, n = 999)
  expect_equal(ncol(result), ncol(se))
})

test_that("select_random_samples is reproducible with set.seed", {
  set.seed(42); r1 <- select_random_samples(se, n = 10)
  set.seed(42); r2 <- select_random_samples(se, n = 10)
  expect_equal(colnames(r1), colnames(r2))
})

# ── prepare_expression_matrix() ───────────────────────────────────────────────

test_that("prepare_expression_matrix errors on non-SummarizedExperiment", {
  expect_error(
    prepare_expression_matrix("not_a_se"),
    "se must be a SummarizedExperiment object"
  )
})

test_that("prepare_expression_matrix errors when counts assay is missing", {
  se_no_counts <- SummarizedExperiment(
    assays = list(tpm = assay(se, "counts"))
  )
  expect_error(
    prepare_expression_matrix(se_no_counts),
    "se must contain a 'counts' assay"
  )
})

test_that("prepare_expression_matrix errors when gene_name column is missing", {
  se_no_gene_name <- se
  rowData(se_no_gene_name)$gene_name <- NULL
  expect_error(
    prepare_expression_matrix(se_no_gene_name),
    "rowData\\(se\\) must contain a 'gene_name' column"
  )
})

test_that("prepare_expression_matrix returns a matrix", {
  result <- prepare_expression_matrix(se)
  expect_true(is.matrix(result))
})

test_that("prepare_expression_matrix returns correct dimensions", {
  result <- prepare_expression_matrix(se)
  expect_equal(nrow(result), nrow(se))
  expect_equal(ncol(result), ncol(se))
})

test_that("prepare_expression_matrix applies log2(x + 1) normalization", {
  result <- prepare_expression_matrix(se)
  raw    <- assay(se, "counts")[1, 1]
  expect_equal(result[1, 1], log2(raw + 1))
})

test_that("prepare_expression_matrix values are non-negative", {
  result <- prepare_expression_matrix(se)
  expect_true(all(result >= 0))
})

test_that("prepare_expression_matrix uses gene symbols as rownames", {
  result <- prepare_expression_matrix(se)
  gene_names <- rowData(se)$gene_name
  expect_true(any(rownames(result) %in% gene_names))
})

test_that("prepare_expression_matrix produces unique rownames", {
  result <- prepare_expression_matrix(se)
  expect_equal(length(rownames(result)), length(unique(rownames(result))))
})
