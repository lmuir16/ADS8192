# Find and subset hub genes

Rank genes by degree and extract top `n` hub genes from network

## Usage

``` r
get_hub_genes(network, modules = NULL, n = 20)
```

## Arguments

- network:

  A network object output by build_network()

- modules:

  A named numeric vector of module assignments

- n:

  Number of top hub genes to return (default: 20)

## Value

A data frame of hub genes ranked by degree (highest to lowest)

## Examples

``` r
data(example_se)
mat <- prepare_expression_matrix(example_se)
cor_mat <- pairwise_corr(mat)
adj_mat <- build_adjacency_mat(cor_mat)
net <- build_network(adj_mat)
gene_mods <- detect_network_modules(net)
#> Modularity: 0.545
#> Modules found: 514
hub_genes <- get_hub_genes(net, modules = gene_mods$modules, n = 20)
```
