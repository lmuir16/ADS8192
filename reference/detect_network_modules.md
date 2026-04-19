# Detect gene modules/communities

Performs community detection on a gene co-expression network using
Louvain or walktrap clustering algorithms and returns module assignments
and modularity statistics.

## Usage

``` r
detect_network_modules(network, community_method = "louvain")
```

## Arguments

- network:

  A network object output by build_network()

- community_method:

  Community detection method: "louvain" or "walktrap" (default:
  "louvain").

## Value

A list containing:

- community:

  An igraph community object

- modules:

  A named numeric vector mapping genes to module IDs

- modularity:

  A numeric modularity score

- module_df:

  A data frame of gene-module assignments with node-level metadata (e.g.
  degree)

## Examples

``` r
data(example_se)
mat <- prepare_expression_matrix(example_se)
cor_mat <- pairwise_corr(mat)
adj_mat <- build_adjacency_mat(cor_mat)
net <- build_network(adj_mat)
gene_mods <- detect_network_modules(net, community_method = "louvain")
#> Modularity: 0.547
#> Modules found: 513
```
