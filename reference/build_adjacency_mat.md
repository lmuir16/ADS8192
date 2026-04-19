# Threshold to adjacency matrix

Converts gene-gene correlation matrix into a binary adjacency matrix
according to threshold for absolute correlation value.

## Usage

``` r
build_adjacency_mat(cor_mat, cor_threshold = 0.7)
```

## Arguments

- cor_mat:

  The gene-gene correlation matrix

- cor_threshold:

  Correlation threshold for defining edges (default: 0.7)

## Value

A binary adjacency matrix with 1 indicating an edge and 0 indicating no
edge

## Examples

``` r
data(example_se)
mat <- prepare_expression_matrix(example_se)
cor_mat <- pairwise_corr(mat)
adj_mat <- build_adjacency_mat(cor_mat, cor_threshold=0.7)
```
