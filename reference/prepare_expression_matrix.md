# Prepare log-normalized expression matrix from a SummarizedExperiment

Extracts the counts assay from a SummarizedExperiment object, applies
log2(x + 1) normalization, and remaps Ensembl row identifiers to HGNC
gene symbols using rowData. Falls back to Ensembl IDs for genes with
missing or empty symbols, and resolves duplicate symbols with
[`make.unique()`](https://rdrr.io/r/base/make.unique.html).

## Usage

``` r
prepare_expression_matrix(se)
```

## Arguments

- se:

  A SummarizedExperiment object with a `counts` assay If
  [`rowData()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  contains a `gene_name` column, Ensembl IDs are remapped to HGNC gene
  symbols. Otherwise rownames are used as-is.

## Value

A log2-normalized genes x samples matrix with gene symbols as row names

## Examples

``` r
library(SummarizedExperiment)
data(example_se)
mat <- prepare_expression_matrix(example_se)
```
