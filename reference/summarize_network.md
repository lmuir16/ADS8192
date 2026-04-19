# Summarize global network properties

Summarize global network properties

## Usage

``` r
summarize_network(network, modules, modularity_score)
```

## Arguments

- network:

  A network object output by build_network()

- modules:

  A named numeric vector of module assignments

- modularity_score:

  A numeric modularity score

## Value

A data frame of global network statistics

## Examples

``` r
data(example_se)
mat <- prepare_expression_matrix(example_se)
cor_mat <- pairwise_corr(mat)
adj_mat <- build_adjacency_mat(cor_mat)
net <- build_network(adj_mat)
gene_mods <- detect_network_modules(net)
#> Modularity: 0.533
#> Modules found: 518
global_net_summary <- summarize_network(net, gene_mods$modules,
    gene_mods$modularity)
#>          metric      value
#> 1       n_genes  2000.0000
#> 2       n_edges 29387.0000
#> 3   mean_degree    29.3900
#> 4 median_degree     7.0000
#> 5  n_hubs_deg20   719.0000
#> 6     n_modules   518.0000
#> 7    modularity     0.5325
```
