# Plot gene-gene correlation heatmap

Plot gene-gene correlation heatmap

## Usage

``` r
plot_corr_map(cor_mat, cor_method = "pearson")
```

## Arguments

- cor_mat:

  A gene-gene correlation matrix

- cor_method:

  Correlation method used to compute matrix (default: "pearson")

## Value

Renders the heatmap to the active graphics device. The ComplexHeatmap
object is returned invisibly for optional reuse.

## Examples

``` r
data(example_se)
mat <- prepare_expression_matrix(example_se)
cor_mat <- pairwise_corr(mat)
ht <- plot_corr_map(cor_mat, cor_method = "pearson")
```
