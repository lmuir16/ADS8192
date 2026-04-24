# tests/testthat/test-export.R
# Tests for net_export()
#
# Shared fixture using a subset of example_se (GTEx data)

data(example_se)
se          <- example_se[1:50, 1:30]
mat         <- prepare_expression_matrix(se)
cor_mat     <- pairwise_corr(mat)
adj_mat     <- build_adjacency_mat(cor_mat, cor_threshold = 0.3)
net         <- build_network(adj_mat)
gene_mods   <- detect_network_modules(net)
net_summary <- summarize_network(net, gene_mods$modules, gene_mods$modularity)

# ── Input validation ──────────────────────────────────────────────────────────

test_that("net_export errors on invalid network input", {
  expect_error(
    net_export(cor_mat, list(graph = "not_igraph"), net_summary,
               gene_mods$modules, gene_mods$module_df),
    "network must contain an igraph object in \\$graph"
  )
})

test_that("net_export errors on missing network$degree", {
  expect_error(
    net_export(cor_mat, list(graph = net$graph), net_summary,
               gene_mods$modules, gene_mods$module_df),
    "network\\$degree is missing"
  )
})

test_that("net_export errors on missing network$betweenness", {
  expect_error(
    net_export(cor_mat,
               list(graph = net$graph, degree = net$degree),
               net_summary, gene_mods$modules, gene_mods$module_df),
    "network\\$betweenness is missing"
  )
})

test_that("net_export warns when n_hubs exceeds network size", {
  expect_warning(
    net_export(cor_mat, net, net_summary, gene_mods$modules,
               gene_mods$module_df, n_hubs = 999),
    "n_hubs"
  )
})

# ── Return value ──────────────────────────────────────────────────────────────

test_that("net_export returns a character vector of file names", {
  result <- net_export(cor_mat, net, net_summary, gene_mods$modules,
                       gene_mods$module_df)
  expect_type(result, "character")
})

test_that("net_export returns all six expected file names", {
  result <- net_export(cor_mat, net, net_summary, gene_mods$modules,
                       gene_mods$module_df)
  expected_files <- c("gene_correlations.tsv", "hub_genes.tsv",
                      "module_assignments.tsv", "module_summary.tsv",
                      "network_summary.tsv", "node_stats.tsv")
  expect_true(all(expected_files %in% result))
})

# ── File contents ─────────────────────────────────────────────────────────────

test_that("net_export writes files that are readable TSVs", {
  net_export(cor_mat, net, net_summary, gene_mods$modules,
             gene_mods$module_df)
  output_dir <- file.path(tempdir(), "network_output")

  for (f in c("gene_correlations.tsv", "hub_genes.tsv",
              "module_assignments.tsv", "module_summary.tsv",
              "network_summary.tsv", "node_stats.tsv")) {
    result <- read.table(file.path(output_dir, f),
                         sep = "\t", header = TRUE)
    expect_s3_class(result, "data.frame")
    expect_gt(nrow(result), 0)
  }
})

test_that("gene_correlations.tsv contains correct columns", {
  net_export(cor_mat, net, net_summary, gene_mods$modules,
             gene_mods$module_df)
  output_dir <- file.path(tempdir(), "network_output")
  result <- read.table(file.path(output_dir, "gene_correlations.tsv"),
                       sep = "\t", header = TRUE)
  expect_named(result, c("gene1", "gene2", "correlation"))
})

test_that("gene_correlations.tsv correlations respect threshold", {
  net_export(cor_mat, net, net_summary, gene_mods$modules,
             gene_mods$module_df, cor_threshold = 0.3)
  output_dir <- file.path(tempdir(), "network_output")
  result <- read.table(file.path(output_dir, "gene_correlations.tsv"),
                       sep = "\t", header = TRUE)
  expect_true(all(abs(result$correlation) >= 0.3))
})

test_that("node_stats.tsv contains correct columns", {
  net_export(cor_mat, net, net_summary, gene_mods$modules,
             gene_mods$module_df)
  output_dir <- file.path(tempdir(), "network_output")
  result <- read.table(file.path(output_dir, "node_stats.tsv"),
                       sep = "\t", header = TRUE)
  expect_named(result, c("gene", "degree", "betweenness", "module"))
})

test_that("node_stats.tsv covers all genes in the network", {
  net_export(cor_mat, net, net_summary, gene_mods$modules,
             gene_mods$module_df)
  output_dir <- file.path(tempdir(), "network_output")
  result <- read.table(file.path(output_dir, "node_stats.tsv"),
                       sep = "\t", header = TRUE)
  expect_equal(nrow(result), length(net$degree))
})

test_that("hub_genes.tsv contains exactly 20 rows by default", {
  net_export(cor_mat, net, net_summary, gene_mods$modules,
             gene_mods$module_df)
  output_dir <- file.path(tempdir(), "network_output")
  result <- read.table(file.path(output_dir, "hub_genes.tsv"),
                       sep = "\t", header = TRUE)
  expect_equal(nrow(result), 20)
})

test_that("hub_genes.tsv respects n_hubs parameter", {
  net_export(cor_mat, net, net_summary, gene_mods$modules,
             gene_mods$module_df, n_hubs = 10)
  output_dir <- file.path(tempdir(), "network_output")
  result <- read.table(file.path(output_dir, "hub_genes.tsv"),
                       sep = "\t", header = TRUE)
  expect_equal(nrow(result), 10)
})

test_that("network_summary.tsv contains all expected metrics", {
  net_export(cor_mat, net, net_summary, gene_mods$modules,
             gene_mods$module_df)
  output_dir <- file.path(tempdir(), "network_output")
  result <- read.table(file.path(output_dir, "network_summary.tsv"),
                       sep = "\t", header = TRUE)
  expected_metrics <- c("n_genes", "n_edges", "mean_degree",
                        "median_degree", "n_hubs_deg20",
                        "n_modules", "modularity")
  expect_true(all(expected_metrics %in% result$metric))
})

test_that("module_summary.tsv contains correct columns", {
  net_export(cor_mat, net, net_summary, gene_mods$modules,
             gene_mods$module_df)
  output_dir <- file.path(tempdir(), "network_output")
  result <- read.table(file.path(output_dir, "module_summary.tsv"),
                       sep = "\t", header = TRUE)
  expect_named(result, c("module", "n_genes", "mean_degree",
                         "mean_betweenness", "top_hub"))
})

test_that("module_summary.tsv gene counts sum to total genes", {
  net_export(cor_mat, net, net_summary, gene_mods$modules,
             gene_mods$module_df)
  output_dir <- file.path(tempdir(), "network_output")
  result <- read.table(file.path(output_dir, "module_summary.tsv"),
                       sep = "\t", header = TRUE)
  expect_equal(sum(result$n_genes), length(net$degree))
})
