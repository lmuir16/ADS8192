# R/plotting.R

# Declaring names as known globals to prevent unnecessary check() warnings
utils::globalVariables(c("x", "y", "xend", "yend", "degree", "module", "is_hub",
                         "gene"))

#' Plot gene-gene correlation heatmap
#'
#' @param cor_mat A gene-gene correlation matrix
#' @param cor_method Correlation method used to compute matrix
#' (default: "pearson")
#'
#' @return Renders the heatmap to the active graphics device.
#'   The ComplexHeatmap object is returned invisibly for optional reuse.
#'
#' @importFrom ComplexHeatmap Heatmap draw
#'
#' @export
#'
#' @examples
#' data(example_se)
#' mat <- prepare_expression_matrix(example_se)
#' cor_mat <- pairwise_corr(mat)
#' ht <- plot_corr_map(cor_mat, cor_method = "pearson")

plot_corr_map <- function(cor_mat, cor_method = "pearson") {

  # Input check
  if (!is.matrix(cor_mat)) {
    stop("cor_mat must be a matrix")
  }

  ht <- Heatmap(cor_mat, name = paste0(tools::toTitleCase(cor_method), " r"),
                show_row_names = FALSE, show_column_names = FALSE,
                column_title = paste0("Gene-Gene Correlation (top ",
                                      nrow(cor_mat), " genes)"))
  draw(ht)
  invisible(ht)
}

#' Plot force-directed network graph of hub genes
#'
#' Plots a force-directed network graph of the top \code{n_top} hub genes
#' colored by module membership. The top 20 hub genes are labeled by gene symbol.
#'
#' @param network A network object output by build_network()
#' @param modules A named integer vector of module assignments
#' @param n_top The number of top hub genes to plot, ranked by degree (default: 80)
#' @param cor_threshold Correlation threshold used to define adjacency (default: 0.70)
#' @param community_method Community detection method used to detect modules
#'   (default: "louvain")
#' @param n_label Number of top hub genes to label by gene symbol (default: 20)
#'
#' @return Renders the network graph to the active graphics device.
#'   The ggplot object is returned invisibly for optional reuse.
#'
#' @importFrom ggplot2 ggplot geom_segment geom_point geom_label aes
#'   scale_size_continuous scale_color_brewer theme_void labs
#' @importFrom grid unit
#' @importFrom igraph induced_subgraph layout_with_fr as_edgelist V "V<-" degree
#' @export
#'
#' @examples
#' data(example_se)
#' mat <- prepare_expression_matrix(example_se)
#' cor_mat <- pairwise_corr(mat)
#' adj_mat <- build_adjacency_mat(cor_mat)
#' net <- build_network(adj_mat)
#' gene_mods <- detect_network_modules(net)
#' plot_network(net, modules = gene_mods$modules, n_top = 80, n_label = 20)

plot_network <- function(network, modules, n_top = 80, cor_threshold = 0.7,
                         community_method = "louvain", n_label = 20) {

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

  # Check that n_top is positive
  if (!is.numeric(n_top) || length(n_top) != 1 || n_top < 1) {
    stop("n_top must be a single positive integer.")
  }

  # Check that n_top <= length(degree)
  if (n_top > length(deg)){
    stop(
      paste0(
        "n_top (", n_top,
        ") exceeds the number of nodes in the graph (",
        length(deg), "). ",
        "Please choose n_top <= ", length(deg), "."
      )
    )
  }

  # Label all plotted hub genes if n_label > n_top (warning)
  if (n_label > n_top) {
    warning(paste0("n_label (", n_label, ") exceeds n_top (", n_top,
                   "). Labeling all ", n_top, " plotted genes."))
    n_label <- n_top
  }

  # Subset to top n hub genes to keep layout readable

  ##  Sort genes highest to lowest by degree and extract top n gene names
  hub_names <- get_hub_genes(network, n = n_top)$gene
  ## Top 20 hub genes
  top_label_names <- head(hub_names, n_label)
  ## Subgraph containing specified vertices and all their edges
  g_sub <- induced_subgraph(g, vids = hub_names)

  # Module membership and degree for the subgraph
  V(g_sub)$module   <- as.character(modules[V(g_sub)$name])
  V(g_sub)$degree   <- degree(g_sub)
  V(g_sub)$is_label <- V(g_sub)$name %in% top_label_names

  # Layout
  layout_fr <- layout_with_fr(g_sub)

  # Build a ggplot-friendly data frame from the igraph layout
  node_df <- data.frame(
    gene        = V(g_sub)$name,   # gene symbol (or Ensembl fallback)
    x           = layout_fr[, 1],
    y           = layout_fr[, 2],
    degree      = V(g_sub)$degree,
    module      = paste0("M", V(g_sub)$module),
    is_hub      = V(g_sub)$is_label,
    stringsAsFactors = FALSE
  )

  edge_list <- as_edgelist(g_sub)
  edge_df <- data.frame(
    x    = layout_fr[match(edge_list[, 1], V(g_sub)$name), 1],
    y    = layout_fr[match(edge_list[, 1], V(g_sub)$name), 2],
    xend = layout_fr[match(edge_list[, 2], V(g_sub)$name), 1],
    yend = layout_fr[match(edge_list[, 2], V(g_sub)$name), 2]
  )

  p_net <- ggplot() +
    geom_segment(data = edge_df,
                 aes(x = x, y = y, xend = xend, yend = yend),
                 color = "grey75", linewidth = 0.25, alpha = 0.5) +
    geom_point(data = node_df,
               aes(x = x, y = y, size = degree, color = module),
               alpha = 0.85) +
    geom_label(data = subset(node_df, is_hub),
               aes(x = x, y = y, label = gene),
               size = 2.5, label.padding = unit(0.15, "lines"),
               linewidth = 0.2, alpha = 0.85) +
    scale_size_continuous(range = c(2, 8), name = "Degree") +
    scale_color_brewer(palette = "Set2", name = "Module") +
    theme_void() +
    labs(title = "Hub Gene Co-expression Network",
         subtitle = paste0("Top ", n_top, " genes by degree | cor threshold = ",
                           cor_threshold, " | ",
                           tools::toTitleCase(community_method), " modules = ",
                           length(unique(modules))))
  print(p_net)
  invisible(p_net)
}
