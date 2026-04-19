# Export network results as TSV files

Writes correlation edges, network summary statistics, hub genes, module
assignments, node statistics, and module summaries to TSV files

## Usage

``` r
net_export(
  cor_mat,
  network,
  net_summary,
  modules,
  module_df,
  cor_threshold = 0.7
)
```

## Arguments

- cor_mat:

  A gene-gene correlation matrix

- network:

  A network object output by build_network()

- net_summary:

  A data frame of global network statistics

- modules:

  A named integer vector mapping genes to module IDs

- module_df:

  A data frame of gene-module assignments with node-level metadata (e.g.
  degree)

- cor_threshold:

  Correlation threshold used to define significant edges (default: 0.7;
  should match the value passed to build_adjacency_mat())

## Value

A character vector of exported file names

## Examples

``` r
data(example_se)
mat <- prepare_expression_matrix(example_se)
cor_mat <- pairwise_corr(mat)
adj_mat <- build_adjacency_mat(cor_mat)
net <- build_network(adj_mat)
gene_mods <- detect_network_modules(net)
#> Modularity: 0.547
#> Modules found: 513
net_summary <- summarize_network(net, gene_mods$modules, gene_mods$modularity)
#>          metric      value
#> 1       n_genes  2000.0000
#> 2       n_edges 29387.0000
#> 3   mean_degree    29.3900
#> 4 median_degree     7.0000
#> 5  n_hubs_deg20   719.0000
#> 6     n_modules   513.0000
#> 7    modularity     0.5472
network_results <- net_export(cor_mat, net, net_summary, gene_mods$modules,
    gene_mods$module_df, cor_threshold = 0.7)
```
