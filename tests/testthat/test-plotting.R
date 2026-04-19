# tests/testthat/test-plotting.R
# Tests for plotting functions
#
# Shared fixture using a subset of example_se (real GTEx data)

data(example_se)
se        <- example_se[1:50, 1:30]
mat       <- prepare_expression_matrix(se)
cor_mat   <- pairwise_corr(mat)
adj_mat   <- build_adjacency_mat(cor_mat, cor_threshold = 0.3)
net       <- build_network(adj_mat)
gene_mods <- detect_network_modules(net)

# ── plot_corr_map() ───────────────────────────────────────────────────────────

test_that("plot_corr_map returns a Heatmap object invisibly", {
  result <- plot_corr_map(cor_mat)
  expect_s4_class(result, "Heatmap")
})

test_that("plot_corr_map accepts cor_method argument", {
  expect_no_error(plot_corr_map(cor_mat, cor_method = "spearman"))
})

test_that("plot_corr_map errors on non-matrix input", {
  expect_error(plot_corr_map("not_a_matrix"))
})

# ── plot_network() ────────────────────────────────────────────────────────────

test_that("plot_network returns a ggplot object invisibly", {
  result <- plot_network(net, modules = gene_mods$modules, n_top = 20)
  expect_s3_class(result, "ggplot")
})

test_that("plot_network errors on invalid network input", {
  expect_error(
    plot_network(list(graph = "not_igraph"), gene_mods$modules),
    "network must contain an igraph object in \\$graph"
  )
})

test_that("plot_network errors on invalid n_top", {
  expect_error(
    plot_network(net, gene_mods$modules, n_top = -1),
    "n_top must be a single positive integer"
  )
})

test_that("plot_network errors when n_top exceeds node count", {
  expect_error(
    plot_network(net, gene_mods$modules, n_top = 999),
    "exceeds the number of nodes"
  )
})

test_that("plot_network errors on missing network$degree", {
  expect_error(
    plot_network(list(graph = net$graph), gene_mods$modules),
    "network\\$degree is missing"
  )
})
