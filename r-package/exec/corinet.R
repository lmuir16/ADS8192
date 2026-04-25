#!/usr/bin/env Rapp
#| name: corinet
#| title: Gene Co-Expression Network Analysis
#| description: Build and analyze gene co-expression networks from bulk RNA-seq data.

suppressPackageStartupMessages({
  library(corinet)
  library(utils)
  library(stats)
  library(ggplot2)
  library(grDevices)
  library(SummarizedExperiment)
})

# Helper to read TSV/CSV (not exported; kept in CLI script)
read_data_file <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "csv") {
    utils::read.csv(path, row.names = 1, check.names = FALSE)
  } else {
    utils::read.table(path, sep = "\t", header = TRUE, row.names = 1,
                      check.names = FALSE)
  }
}

message("

 ▄▄▄▄  ▄▄▄  ▄▄▄▄  ▄▄ ▄▄  ▄▄ ▄▄▄▄▄ ▄▄▄▄▄▄
██▀▀▀ ██▀██ ██▄█▄ ██ ███▄██ ██▄▄    ██
▀████ ▀███▀ ██ ██ ██ ██ ▀██ ██▄▄▄   ██






                                         ")

switch(
  "",

  #| title: Run gene network analysis pipeline
  #| description: Build a gene co-expression network from a counts matrix, export results.
  network = {

    #| description: Path to counts matrix (TSV/CSV, genes x samples)
    #| short: c
    counts <- ""

    #| description: Output directory
    #| short: o
    output <- ""

    #| description: Number of top variable genes to use
    #| short: g
    n_top <- 500L

    #| description: Correlation method to use (pearson or spearman)
    #| short: m
    cor_method <- "pearson"

    #| description: Correlation threshold for defining edges
    #| short: t
    cor_threshold <- 0.7

    #| description: Community detection method (louvain or walktrap)
    #| short: d
    community_method <- "louvain"

    #| description: Number of hub genes to export
    #| short: h
    n_hubs <- 20L

    # Validation

    # required arguments (no default)
    if (counts == "" || output == "") {
      stop("--counts and --output are required", call. = FALSE)
    }

    # counts file exists
    if (!file.exists(counts)) {
      stop("File not found: ", counts, call. = FALSE)
    }

    # create output directory if it does not already exist
    if (!dir.exists(output)) dir.create(output, recursive = TRUE)

    # Read counts matrix
    counts_df <- read_data_file(counts)

    # Build SummarizedExperiment
    se <- SummarizedExperiment(
      assays = list(counts = as.matrix(counts_df))
    )

    # Run core network analysis pipeline
    se_top    <- top_variable_features(se, n = n_top)
    mat       <- prepare_expression_matrix(se_top)
    cor_mat   <- pairwise_corr(mat, cor_method = cor_method)
    adj_mat   <- build_adjacency_mat(cor_mat, cor_threshold = cor_threshold)
    net       <- build_network(adj_mat)
    gene_mods <- detect_network_modules(net,
                                        community_method = community_method)
    net_summary <- summarize_network(net, gene_mods$modules,
                                     gene_mods$modularity)

    # Export results
    exported <- net_export(
      cor_mat       = cor_mat,
      cor_threshold = cor_threshold,
      network       = net,
      net_summary   = net_summary,
      modules       = gene_mods$modules,
      module_df     = gene_mods$module_df,
      n_hubs        = n_hubs,
      output_dir    = output
    )

    message("Huzzah! Analysis complete. ", length(exported),
            " files written to ", output)

  },

  #| title: Plot correlation heatmap
  #| description: Generate a gene-gene correlation heatmap from a counts matrix.
  heatmap = {

    #| description: Path to counts matrix (TSV/CSV, genes x samples)
    #| short: c
    counts <- ""

    #| description: Output directory
    #| short: o
    output <- ""

    #| description: Number of top variable genes to use
    #| short: g
    n_top <- 500L

    #| description: Correlation method to use (pearson or spearman)
    #| short: m
    cor_method <- "pearson"

    # Validation

    # required arguments (no default)
    if (counts == "" || output == "") {
      stop("--counts and --output are required", call. = FALSE)
    }

    # counts file exists
    if (!file.exists(counts)) {
      stop("File not found: ", counts, call. = FALSE)
    }

    # create output directory if it does not already exist
    if (!dir.exists(output)) dir.create(output, recursive = TRUE)

    # Read counts matrix
    counts_df <- read_data_file(counts)

    # Build SummarizedExperiment
    se <- SummarizedExperiment(
      assays = list(counts = as.matrix(counts_df))
    )

    # Prepare correlation matrix
    se_top <- top_variable_features(se, n = n_top)
    mat <- prepare_expression_matrix(se_top)
    cor_mat <- pairwise_corr(mat, cor_method = cor_method)

    # Plot and save in output directory
    plot_file <- file.path(output, "correlation_heatmap.png")
    grDevices::png(plot_file, width = 10, height = 9, units = "in", res = 300)
    plot_corr_map(cor_mat, cor_method = cor_method)
    dev.off() # closes PNG graphics device

    message("Saved: ", plot_file)
    message("All done.")

  },

  #| title: Plot network graph of hub genes
  #| description: Generate a force-directed network graph from a counts matrix.
  plot_network = {

    #| description: Path to counts matrix (TSV/CSV, genes x samples)
    #| short: c
    counts <- ""

    #| description: Output directory
    #| short: o
    output <- ""

    #| description: Number of top variable genes to use
    #| short: g
    n_top <- 500L

    #| description: Correlation method to use (pearson or spearman)
    #| short: m
    cor_method <- "pearson"

    #| description: Correlation threshold for defining edges
    #| short: t
    cor_threshold <- 0.7

    #| description: Community detection method (louvain or walktrap)
    #| short: d
    community_method <- "louvain"

    #| description: Number of top hub genes to plot
    #| short: k
    n_plot <- 80L

    #| description: Number of top hub genes to label by gene symbol
    #| short: l
    n_label <- 20L

    # Validation

    # required arguments (no default)
    if (counts == "" || output == "") {
      stop("--counts and --output are required", call. = FALSE)
    }

    # counts file exists
    if (!file.exists(counts)) {
      stop("File not found: ", counts, call. = FALSE)
    }

    # create output directory if it does not already exist
    if (!dir.exists(output)) dir.create(output, recursive = TRUE)

    # Read counts matrix
    counts_df <- read_data_file(counts)

    # Build SummarizedExperiment
    se <- SummarizedExperiment(
      assays = list(counts = as.matrix(counts_df))
    )

    # Run data preparation and network analysis
    se_top    <- top_variable_features(se, n = n_top)
    mat       <- prepare_expression_matrix(se_top)
    cor_mat   <- pairwise_corr(mat, cor_method = cor_method)
    adj_mat   <- build_adjacency_mat(cor_mat, cor_threshold = cor_threshold)
    net       <- build_network(adj_mat)
    gene_mods <- detect_network_modules(net,
                                        community_method = community_method)

    # Plot and save in output directory
    plot_file <- file.path(output, "network_plot.png")
    grDevices::png(plot_file, width = 10, height = 9, units = "in", res = 300)
    plot_network(net,
                 modules = gene_mods$modules,
                 n_top = n_plot,
                 cor_threshold = cor_threshold,
                 community_method = community_method)
    dev.off() # closes PNG graphics device

    message("Saved: ", plot_file)
    message("All done.")

  }

)
