# R/export.R

#' Export network results as TSV files
#'
#' Writes correlation edges, network summary statistics, hub genes,
#' module assignments, node statistics, and module summaries to TSV files
#'
#' @param cor_mat A gene-gene correlation matrix
#' @param network A network object output by build_network()
#' @param net_summary A data frame of global network statistics
#' @param modules A named integer vector mapping genes to module IDs
#' @param module_df A data frame of gene-module assignments
#'    with node-level metadata (e.g. degree)
#' @param cor_threshold Correlation threshold used to define significant edges
#'    (default: 0.7; should match the value passed to build_adjacency_mat())
#' @param n_hubs Number of top hub genes to export to hub_genes.tsv (default: 20)
#' @param output_dir Directory to write TSV files to
#'    (default: a "network_output" folder in the session temp directory)
#'
#' @return A character vector of exported file names
#'
#' @importFrom utils write.table
#' @export
#'
#' @examples
#' data(example_se)
#' mat <- prepare_expression_matrix(example_se)
#' cor_mat <- pairwise_corr(mat)
#' adj_mat <- build_adjacency_mat(cor_mat)
#' net <- build_network(adj_mat)
#' gene_mods <- detect_network_modules(net)
#' net_summary <- summarize_network(net, gene_mods$modules, gene_mods$modularity)
#' network_results <- net_export(cor_mat, net, net_summary, gene_mods$modules,
#'     gene_mods$module_df, cor_threshold = 0.7, n_hubs = 10)

net_export <- function(cor_mat, network, net_summary, modules,
                       module_df, cor_threshold = 0.7, n_hubs = 20,
                       output_dir = file.path(tempdir(), "network_output")) {

  # Input checks for network
  if (!is.list(network) || !inherits(network$graph, "igraph")) {
    stop("network must contain an igraph object in $graph")
  }

  if (is.null(network$degree)) {
    stop(paste0("network$degree is missing. ",
                "Ensure network was created by build_network()."))
  }

  if (is.null(network$betweenness)) {
    stop(paste0("network$betweenness is missing. ",
                "Ensure network was created by build_network()."))
  }

  # Return all available genes if n_hub > n_genes (warning)
  n_genes <- length(network$degree)
  if (n_hubs > n_genes) {
    warning(paste0("n_hubs (", n_hubs, ") exceeds the number of genes in the ",
                   "network (", n_genes, "). Exporting all ", n_genes, " genes."))
    n_hubs <- n_genes
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Write top correlations (flattened, filtered to |r| >= threshold)
  # rownames are gene symbols after the remap above
  upper_idx <- which(upper.tri(cor_mat) & abs(cor_mat) >= cor_threshold,
                     arr.ind = TRUE)
  cor_edges <- data.frame(
    gene1       = rownames(cor_mat)[upper_idx[, 1]],
    gene2       = colnames(cor_mat)[upper_idx[, 2]],
    correlation = cor_mat[upper_idx]
  )

  g <- network$graph
  deg <- network$degree
  betw <- network$betweenness
  hub_genes <- get_hub_genes(network, modules, n = n_hubs)

  write.table(cor_edges, file.path(output_dir, "gene_correlations.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(net_summary, file.path(output_dir, "network_summary.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(hub_genes, file.path(output_dir, "hub_genes.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(module_df, file.path(output_dir, "module_assignments.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  # All-node statistics (degree + betweenness + module for every gene)
  node_stats <- data.frame(
    gene        = names(deg),
    degree      = deg,
    betweenness = betw,
    module      = paste0("M", modules[names(deg)])
  )
  node_stats <- node_stats[order(node_stats$degree, decreasing = TRUE), ]
  write.table(node_stats, file.path(output_dir, "node_stats.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  # Per-module summary: size, mean degree, mean betweenness
  mod_summary <- do.call(rbind, lapply(unique(node_stats$module), function(m) {
    sub <- node_stats[node_stats$module == m, ]
    data.frame(
      module = m,
      n_genes = nrow(sub),
      mean_degree = round(mean(sub$degree), 2),
      mean_betweenness = round(mean(sub$betweenness), 2),
      top_hub = sub$gene[which.max(sub$degree)]
    )
  }))
  mod_summary <- mod_summary[order(mod_summary$n_genes, decreasing = TRUE), ]
  write.table(mod_summary, file.path(output_dir, "module_summary.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  list.files(output_dir)
}
