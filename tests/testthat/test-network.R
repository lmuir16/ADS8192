# tests/testthat/test-network.R
# Tests for network analysis functions
#
# Shared fixture using a subset of example_se (GTEx data) for computational
# efficiency -- 50 genes x 30 samples

data(example_se)
se        <- example_se[1:50, 1:30]
mat       <- prepare_expression_matrix(se)
cor_mat   <- pairwise_corr(mat)
adj_mat   <- build_adjacency_mat(cor_mat, cor_threshold = 0.3)
net       <- build_network(adj_mat)
gene_mods <- detect_network_modules(net)

# ── pairwise_corr() ───────────────────────────────────────────────────────────

test_that("pairwise_corr errors when passed non-matrix", {
  expect_error(
    pairwise_corr("not_a_matrix"),
    "mat must be a matrix"
  )
})

test_that("pairwise_corr returns a square symmetric matrix", {
  expect_true(is.matrix(cor_mat))
  expect_equal(nrow(cor_mat), ncol(cor_mat))
  expect_equal(cor_mat, t(cor_mat))
})

test_that("pairwise_corr values are bounded [-1, 1]", {
  expect_gte(min(cor_mat, na.rm = TRUE), -1)
  expect_lte(max(cor_mat, na.rm = TRUE),  1)
})

test_that("pairwise_corr diagonal is 1", {
  expect_true(all(diag(cor_mat) == 1))
})

test_that("pairwise_corr respects cor_method argument", {
  spearman <- pairwise_corr(mat, cor_method = "spearman")
  expect_false(identical(cor_mat, spearman))
})

# ── build_adjacency_mat() ─────────────────────────────────────────────────────

test_that("build_adjacency_mat returns a binary symmetric matrix", {
  expect_true(is.matrix(adj_mat))
  expect_equal(dim(adj_mat), dim(cor_mat))
  expect_true(all(adj_mat %in% 0:1))
  expect_equal(adj_mat, t(adj_mat))
})

test_that("build_adjacency_mat diagonal is zero", {
  expect_true(all(diag(adj_mat) == 0))
})

test_that("build_adjacency_mat higher threshold produces fewer edges", {
  adj_low  <- build_adjacency_mat(cor_mat, cor_threshold = 0.2)
  adj_high <- build_adjacency_mat(cor_mat, cor_threshold = 0.8)
  expect_gt(sum(adj_low), sum(adj_high))
})

test_that("build_adjacency_mat errors on non-square matrix", {
  expect_error(
    build_adjacency_mat(matrix(1:6, nrow = 2, ncol = 3)),
    "cor_mat must be square"
  )
})

test_that("build_adjacency_mat warns when no edges detected", {
  expect_warning(
    build_adjacency_mat(cor_mat, cor_threshold = 1.1),
    "No edges detected"
  )
})

test_that("build_adjacency_mat warns when network density is too high", {
  expect_warning(
    build_adjacency_mat(cor_mat, cor_threshold = 0.0),
    "Network density is"
  )
})

# ── build_network() ───────────────────────────────────────────────────────────

test_that("build_network returns correct structure", {
  expect_type(net, "list")
  expect_named(net, c("graph", "degree", "betweenness"))
  expect_s3_class(net$graph, "igraph")
  expect_true(is.numeric(net$degree))
  expect_true(is.numeric(net$betweenness))
})

test_that("build_network degree length matches number of genes", {
  expect_equal(length(net$degree), nrow(adj_mat))
})

test_that("build_network errors on non-square matrix", {
  expect_error(
    build_network(matrix(1:6, nrow = 2)),
    "adj_mat must be a square adjacency matrix"
  )
})

# ── detect_network_modules() ──────────────────────────────────────────────────

test_that("detect_network_modules returns correct structure", {
  expect_type(gene_mods, "list")
  expect_named(gene_mods, c("community", "modules", "modularity", "module_df"))
  expect_s3_class(gene_mods$community, "communities")
  expect_true(is.numeric(gene_mods$modules))
  expect_false(is.null(names(gene_mods$modules)))
  expect_true(is.numeric(gene_mods$modularity))
  expect_s3_class(gene_mods$module_df, "data.frame")
})

test_that("detect_network_modules assigns every gene to a module", {
  expect_equal(length(gene_mods$modules), length(net$degree))
})

test_that("detect_network_modules modularity is between 0 and 1", {
  expect_gte(gene_mods$modularity, 0)
  expect_lte(gene_mods$modularity, 1)
})

test_that("detect_network_modules works with community_method = 'walktrap'", {
  result <- detect_network_modules(net, community_method = "walktrap")
  expect_type(result, "list")
  expect_named(result, c("community", "modules", "modularity", "module_df"))
  expect_s3_class(result$community, "communities")
})

test_that("detect_network_modules errors when passed non-igraph object", {
  expect_error(
    detect_network_modules(list(graph = "not_igraph")),
    regexp = "network must contain an igraph object in \\$graph"
  )
})

test_that("detect_network_modules errors on invalid community_method", {
  expect_error(
    detect_network_modules(net, community_method = "invalid"),
    "community_method must be 'louvain' or 'walktrap'"
  )
})

# ── summarize_network() ───────────────────────────────────────────────────────

test_that("summarize_network returns correct structure", {
  result <- summarize_network(net, gene_mods$modules, gene_mods$modularity)
  expect_s3_class(result, "data.frame")
  expect_true(all(c("n_genes", "n_edges", "mean_degree",
                    "median_degree", "n_hubs_deg20",
                    "n_modules", "modularity") %in% result$metric))
  expect_true(is.numeric(result$value))
})

test_that("summarize_network n_genes matches node count", {
  result <- summarize_network(net, gene_mods$modules, gene_mods$modularity)
  expect_equal(
    result$value[result$metric == "n_genes"],
    length(net$degree)
  )
})

test_that("summarize_network errors on missing network$degree", {
  expect_error(
    summarize_network(list(graph = net$graph), gene_mods$modules,
                      gene_mods$modularity),
    "network\\$degree is missing"
  )
})

# ── get_hub_genes() ───────────────────────────────────────────────────────────

test_that("get_hub_genes returns correct structure with modules", {
  result <- get_hub_genes(net, modules = gene_mods$modules, n = 10)
  expect_s3_class(result, "data.frame")
  expect_named(result, c("gene", "degree", "betweenness", "module"))
  expect_true(is.character(result$gene))
  expect_true(is.character(result$module))
  expect_true(is.numeric(result$degree))
  expect_true(is.numeric(result$betweenness))
  expect_equal(nrow(result), 10)
})

test_that("get_hub_genes returns correct structure without modules", {
  result <- get_hub_genes(net, modules = NULL, n = 10)
  expect_s3_class(result, "data.frame")
  expect_named(result, c("gene", "degree", "betweenness"))
  expect_equal(nrow(result), 10)
})

test_that("get_hub_genes returns genes sorted by degree descending", {
  result <- get_hub_genes(net, n = 20)
  expect_true(all(diff(result$degree) <= 0))
})

test_that("get_hub_genes respects n parameter", {
  expect_equal(nrow(get_hub_genes(net, n = 5)),  5)
  expect_equal(nrow(get_hub_genes(net, n = 30)), 30)
})

test_that("get_hub_genes errors on invalid n", {
  expect_error(
    get_hub_genes(net, n = -1),
    "n must be a single positive integer"
  )
})

test_that("get_hub_genes errors on missing network$betweenness", {
  expect_error(
    get_hub_genes(list(graph = net$graph, degree = net$degree), n = 10),
    "network\\$betweenness is missing"
  )
})

# Huzzah! Tests passed!
