# corinet: Gene Co-expression Network Analysis for Bulk RNA-seq Data

Provides functions to build and analyze gene co-expression networks from
bulk RNA-seq data. Includes tools for selecting the most variable genes,
computing pairwise gene-gene correlation matrices, and thresholding
correlations into binary adjacency matrices with per-gene (node-level)
statistics. Supports community detection using Louvain or walktrap
algorithms with per-module summaries and global network statistics.
Provides visualizations as correlation heatmaps or force-directed
layouts colored by module with hub gene labels. Exports results to TSV
files.

## See also

Key functions:

- [`prepare_expression_matrix()`](prepare_expression_matrix.md) —
  prepare log-normalized matrix from a SummarizedExperiment

- [`top_variable_features()`](top_variable_features.md) — select most
  variable genes

- [`select_random_samples()`](select_random_samples.md) — subset to
  random samples

- [`pairwise_corr()`](pairwise_corr.md) — compute gene-gene correlation
  matrix

- [`build_adjacency_mat()`](build_adjacency_mat.md) — threshold
  correlations to adjacency matrix

- [`build_network()`](build_network.md) — construct igraph network and
  compute node statistics

- [`detect_network_modules()`](detect_network_modules.md) — detect gene
  modules via Louvain or walktrap

- [`summarize_network()`](summarize_network.md) — summarize global
  network properties

- [`get_hub_genes()`](get_hub_genes.md) — extract top hub genes by
  degree

- [`plot_corr_map()`](plot_corr_map.md) — plot correlation heatmap

- [`plot_network()`](plot_network.md) — plot force-directed network
  graph

- [`net_export()`](net_export.md) — export results to TSV files

## Author

**Maintainer**: Liane Muir <liane.muir@stjude.org>
([ORCID](https://orcid.org/0009-0002-5250-0464))
