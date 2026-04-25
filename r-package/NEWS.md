# corinet 0.1.0

## Initial release

* Full gene co-expression network analysis pipeline from bulk RNA-seq data
* Data preparation functions: `top_variable_features()`, `select_random_samples()`, 
  `prepare_expression_matrix()`
* Network construction: `pairwise_corr()`, `build_adjacency_mat()`, `build_network()`
* Community detection: `detect_network_modules()` with Louvain and walktrap support
* Network summarization: `summarize_network()`, `get_hub_genes()`
* Visualization: `plot_corr_map()`, `plot_network()` with module coloring and hub gene labels
* Export: `net_export()` writing six TSV files
* CLI via Rapp with `network`, `heatmap`, and `plot-network` subcommands
* GTEx skeletal muscle example dataset (`example_se`) via recount3
* pkgdown site with getting started vignette
