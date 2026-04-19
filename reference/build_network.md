# Build network from adjacency matrix and compute node-level statistics

Constructs an igraph network object and per-gene (node-level) centrality
measures (degree and betweenness) for downstream analysis.

## Usage

``` r
build_network(adj_mat)
```

## Arguments

- adj_mat:

  A square binary adjacency matrix (genes x genes)

## Value

A list containing:

- graph:

  An igraph object constructed from the adjacency matrix

- degree:

  A named numeric vector of node degrees

- betweenness:

  A named numeric vector of betweenness centralities

## Examples

``` r
data(example_se)
mat <- prepare_expression_matrix(example_se)
cor_mat <- pairwise_corr(mat)
adj_mat <- build_adjacency_mat(cor_mat)
net <- build_network(adj_mat)
```
