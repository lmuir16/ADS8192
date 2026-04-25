# Changelog

## corinet 0.1.0

### Initial release

- Full gene co-expression network analysis pipeline from bulk RNA-seq
  data
- Data preparation functions:
  [`top_variable_features()`](../reference/top_variable_features.md),
  [`select_random_samples()`](../reference/select_random_samples.md),
  [`prepare_expression_matrix()`](../reference/prepare_expression_matrix.md)
- Network construction:
  [`pairwise_corr()`](../reference/pairwise_corr.md),
  [`build_adjacency_mat()`](../reference/build_adjacency_mat.md),
  [`build_network()`](../reference/build_network.md)
- Community detection:
  [`detect_network_modules()`](../reference/detect_network_modules.md)
  with Louvain and walktrap support
- Network summarization:
  [`summarize_network()`](../reference/summarize_network.md),
  [`get_hub_genes()`](../reference/get_hub_genes.md)
- Visualization: [`plot_corr_map()`](../reference/plot_corr_map.md),
  [`plot_network()`](../reference/plot_network.md) with module coloring
  and hub gene labels
- Export: [`net_export()`](../reference/net_export.md) writing six TSV
  files
- CLI via Rapp with `network`, `heatmap`, and `plot-network` subcommands
- GTEx skeletal muscle example dataset (`example_se`) via recount3
- pkgdown site with getting started vignette
