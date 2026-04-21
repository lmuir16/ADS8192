# Pairwise gene-gene correlation matrix

Computes pairwise correlations between genes (rows) across samples
(columns) using the specified correlation method.

## Usage

``` r
pairwise_corr(mat, cor_method = "pearson")
```

## Arguments

- mat:

  A genes x samples matrix (numeric gene expression values)

- cor_method:

  Correlation method to use (default: "pearson")

## Value

A pairwise gene-gene correlation matrix

## Examples

``` r
data(example_se)
mat <- prepare_expression_matrix(example_se)
cor_mat <- pairwise_corr(mat)
```
