# Package index

## Data Preparation

Functions for loading and preparing expression data

- [`top_variable_features()`](top_variable_features.md) : Select top
  variable features
- [`select_random_samples()`](select_random_samples.md) : Select set of
  random samples
- [`prepare_expression_matrix()`](prepare_expression_matrix.md) :
  Prepare log-normalized expression matrix from a SummarizedExperiment

## Network Construction

Functions for building and analyzing co-expression networks

- [`pairwise_corr()`](pairwise_corr.md) : Pairwise gene-gene correlation
  matrix
- [`build_adjacency_mat()`](build_adjacency_mat.md) : Threshold to
  adjacency matrix
- [`build_network()`](build_network.md) : Build network from adjacency
  matrix and compute node-level statistics
- [`detect_network_modules()`](detect_network_modules.md) : Detect gene
  modules/communities
- [`summarize_network()`](summarize_network.md) : Summarize global
  network properties
- [`get_hub_genes()`](get_hub_genes.md) : Find and subset hub genes

## Visualization

Functions for plotting network results

- [`plot_corr_map()`](plot_corr_map.md) : Plot gene-gene correlation
  heatmap
- [`plot_network()`](plot_network.md) : Plot force-directed network
  graph of hub genes

## Export

Functions for exporting results to TSV files

- [`net_export()`](net_export.md) : Export network results as TSV files

## Data

Example datasets

- [`example_se`](example_se.md) : GTEx Skeletal Muscle
  SummarizedExperiment for testing
