# R/network.R

#' Pairwise correlation matrix
#'
#' Computes pairwise correlations between genes (rows) across samples (columns)
#' using the specified correlation method.
#'
#' @param mat A genes x samples matrix (numeric gene expression values)
#' @param cor_method Correlation method to use (default: "pearson")
#'
#' @return A pairwise gene-gene correlation matrix
#'
#' @importFrom stats cor
#' @export
#'
#' @examples
#' data(example_se)
#' mat <- prepare_expression_matrix(example_se)
#' cor_mat <- pairwise_corr(mat)

pairwise_corr <- function(mat, cor_method = "pearson") {

  if (!is.matrix(mat)) {
    stop("mat must be a matrix")
  }

  cor_mat <- cor(t(mat), method = cor_method)

  return(cor_mat)
}

#' Threshold to adjacency matrix
#'
#' Converts gene-gene correlation matrix into a binary adjacency matrix
#' according to threshold for absolute correlation value.
#'
#' @param cor_mat The gene-gene correlation matrix
#' @param cor_threshold Correlation threshold for defining edges (default: 0.7)
#'
#' @return A binary adjacency matrix with 1 indicating an edge
#'   and 0 indicating no edge
#'
#' @export
#'
#' @examples
#' data(example_se)
#' mat <- prepare_expression_matrix(example_se)
#' cor_mat <- pairwise_corr(mat)
#' adj_mat <- build_adjacency_mat(cor_mat, cor_threshold=0.7)

build_adjacency_mat <- function(cor_mat, cor_threshold = 0.7) {
  adj_mat <- (abs(cor_mat) >= cor_threshold) * 1
  diag(adj_mat) <- 0

  return(adj_mat)
}

#' Build network from adjacency matrix and compute node-level statistics
#'
#' Constructs an igraph network object and per-gene (node-level) centrality measures
#' (degree and betweenness) for downstream analysis.
#'
#' @param adj_mat A square binary adjacency matrix (genes x genes)
#' @param graph_mode How igraph should interpret the adjacency matrix
#'    (default: "undirected")
#'
#' @return A list containing:
#' \describe{
#'   \item{graph}{An igraph object constructed from the adjacency matrix}
#'   \item{degree}{A named numeric vector of node degrees}
#'   \item{betweenness}{A named numeric vector of betweenness centralities}
#' }
#'
#' @importFrom igraph graph_from_adjacency_matrix degree betweenness
#' @export
#'
#' @examples
#' data(example_se)
#' mat <- prepare_expression_matrix(example_se)
#' cor_mat <- pairwise_corr(mat)
#' adj_mat <- build_adjacency_mat(cor_mat)
#' net <- build_network(adj_mat)

build_network <- function(adj_mat, graph_mode = "undirected") {

  # Input check for adjacency matrix
  if (!is.matrix(adj_mat) || nrow(adj_mat) != ncol(adj_mat)) {
    stop("adj_mat must be a square adjacency matrix")
  }

  # Construct igraph network
  g <- graph_from_adjacency_matrix(adj_mat, mode = graph_mode)

  list(
    graph = g,
    degree = degree(g),
    betweenness = betweenness(g)
  )
}

#' Detect gene modules/communities
#'
#' Performs community detection on a gene co-expression network using
#' Louvain or walktrap clustering algorithms and
#' returns module assignments and modularity statistics.
#'
#' @param network A network object output by build_network()
#' @param community_method Community detection method: "louvain" or "walktrap"
#'    (default: "louvain")
#'
#' @return A list containing:
#' \describe{
#'   \item{community}{An igraph community object}
#'   \item{modules}{A named integer vector mapping genes to module IDs}
#'   \item{modularity}{A numeric modularity score}
#'   \item{module_df}{A data frame of gene-module assignments with node-level metadata (e.g. degree)}
#' }
#'
#' @importFrom igraph cluster_louvain cluster_walktrap membership modularity
#' @export
#'
#' @examples
#' data(example_se)
#' mat <- prepare_expression_matrix(example_se)
#' cor_mat <- pairwise_corr(mat)
#' adj_mat <- build_adjacency_mat(cor_mat)
#' net <- build_network(adj_mat)
#' gene_mods <- detect_network_modules(net, community_method = "louvain")

detect_network_modules <- function(network, community_method = "louvain") {

  # Input checks

  # For network and degree
  if (!is.list(network) || !inherits(network$graph, "igraph")) {
    stop("network must contain an igraph object in $graph")
  }

  if (is.null(network$degree)) {
    stop(paste0("network$degree is missing. ",
                "Ensure network was created by build_network()."))
  }

  ## for community_method
  if (!community_method %in% c("louvain", "walktrap")) {
    stop("community_method must be 'louvain' or 'walktrap'")
  }

  # Extract graph and degrees from network object
  g <- network$graph
  deg <- network$degree

  # Detect communities
  if (community_method == "louvain"){
    comm <- cluster_louvain(g)
  } else if (community_method == "walktrap"){
    comm <- cluster_walktrap(g)
  }

  # Compute and summarize module membership and modularity
  module_membership <- membership(comm)
  modularity_score <- modularity(comm)

  message("Modularity: ", round(modularity_score, 3))
  message("Modules found: ", length(unique(module_membership)))

  # Alignment safety check for gene names
  if (!identical(sort(names(module_membership)), sort(names(deg)))) {
    warning("Gene names differ between module membership and degree; using intersection for safe alignment")
  }

  # Restrict alignment to shared genes only
  common_genes <- intersect(names(module_membership), names(deg))

  # Build gene-module assignment table
  module_assignments <- data.frame(
    gene   = common_genes,
    module = paste0("M", module_membership[common_genes]),
    degree = deg[common_genes]
  )

  list(
    community = comm,
    modules = module_membership,
    modularity = modularity_score,
    module_df = module_assignments
  )
}

#' Summarize global network properties
#'
#' @param network A network object output by build_network()
#' @param modules A named integer vector of module assignments
#' @param modularity_score A numeric modularity score
#'
#'
#' @return A data frame of global network statistics
#'
#' @importFrom igraph vcount ecount
#' @importFrom stats median
#' @export
#'
#' @examples
#' data(example_se)
#' mat <- prepare_expression_matrix(example_se)
#' cor_mat <- pairwise_corr(mat)
#' adj_mat <- build_adjacency_mat(cor_mat)
#' net <- build_network(adj_mat)
#' gene_mods <- detect_network_modules(net)
#' global_net_summary <- summarize_network(net, gene_mods$modules,
#'     gene_mods$modularity)

summarize_network <- function(network, modules, modularity_score) {

  # Input check for network
  if (!is.list(network) || !inherits(network$graph, "igraph")) {
    stop("network must contain an igraph object in $graph")
  }

  if (is.null(network$degree)) {
    stop(paste0("network$degree is missing. ",
                "Ensure network was created by build_network()."))
  }

  g <- network$graph
  deg <- network$degree

  net_summary <- data.frame(
    metric = c("n_genes", "n_edges", "mean_degree",
               "median_degree", "n_hubs_deg20",
               "n_modules", "modularity"),
    value  = c(vcount(g), ecount(g), round(mean(deg), 2),
               median(deg), sum(deg >= 20),
               length(unique(modules)), round(modularity_score, 4))
  )
  print(net_summary)
  return(net_summary)
}

#' Find and subset hub genes
#'
#' Rank genes by degree and extract top \code{n} hub genes from network
#'
#' @param network A network object output by build_network()
#' @param modules A named integer vector of module assignments
#' @param n Number of top hub genes to return (default: 20)
#'
#' @return A data frame of hub genes ranked by degree (highest to lowest)
#'
#' @importFrom utils head
#' @export
#'
#' @examples
#' data(example_se)
#' mat <- prepare_expression_matrix(example_se)
#' cor_mat <- pairwise_corr(mat)
#' adj_mat <- build_adjacency_mat(cor_mat)
#' net <- build_network(adj_mat)
#' gene_mods <- detect_network_modules(net)
#' hub_genes <- get_hub_genes(net, modules = gene_mods$modules, n = 20)

get_hub_genes <- function(network, modules = NULL, n = 20) {

  # Input checks
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

  if (!is.numeric(n) || length(n) != 1 || n < 1) {
    stop("n must be a single positive integer")
  }

  g <- network$graph
  deg <- network$degree
  betw <- network$betweenness

  # Build data frame with degree and betweenness per gene
  hub_df <- data.frame(
    gene = names(deg),
    degree = deg,
    betweenness = betw,
    stringsAsFactors = FALSE
  )

  # Adds modules assignments column to hub_df if provided
  if (!is.null(modules)) {
    module_vec <- modules[names(deg)]

    if (anyNA(module_vec)) {
      warning("Some genes are missing module assignments; assigning NA.")
    }

    hub_df$module <- ifelse(is.na(module_vec), NA, paste0("M", module_vec))
  }

  # Rank genes by degree
  hub_df <- hub_df[order(hub_df$degree, decreasing = TRUE), ]

  # Subset of top n hub genes
  head(hub_df, n)
}
