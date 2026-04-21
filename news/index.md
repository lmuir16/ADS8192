# Changelog

## ADS8192 0.0.0.9000

### Initial release

- Full gene co-expression network analysis pipeline from bulk RNA-seq
  data
- Data preparation functions:
  [`top_variable_features()`](https://lmuir16.github.io/ADS8192/reference/top_variable_features.md),
  [`select_random_samples()`](https://lmuir16.github.io/ADS8192/reference/select_random_samples.md),
  [`prepare_expression_matrix()`](https://lmuir16.github.io/ADS8192/reference/prepare_expression_matrix.md)
- Network construction:
  [`pairwise_corr()`](https://lmuir16.github.io/ADS8192/reference/pairwise_corr.md),
  [`build_adjacency_mat()`](https://lmuir16.github.io/ADS8192/reference/build_adjacency_mat.md),
  [`build_network()`](https://lmuir16.github.io/ADS8192/reference/build_network.md)
- Community detection:
  [`detect_network_modules()`](https://lmuir16.github.io/ADS8192/reference/detect_network_modules.md)
  with Louvain and walktrap support
- Network summarization:
  [`summarize_network()`](https://lmuir16.github.io/ADS8192/reference/summarize_network.md),
  [`get_hub_genes()`](https://lmuir16.github.io/ADS8192/reference/get_hub_genes.md)
- Visualization:
  [`plot_corr_map()`](https://lmuir16.github.io/ADS8192/reference/plot_corr_map.md),
  [`plot_network()`](https://lmuir16.github.io/ADS8192/reference/plot_network.md)
  with module coloring and hub gene labels
- Export:
  [`net_export()`](https://lmuir16.github.io/ADS8192/reference/net_export.md)
  writing six TSV files
- CLI via Rapp with `network`, `heatmap`, and `plot-network` subcommands
- GTEx skeletal muscle example dataset (`example_se`) via recount3
- pkgdown site with getting started vignette
