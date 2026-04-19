# Select top variable features

Select top variable features

## Usage

``` r
top_variable_features(se, n = 500, assay_name = "counts")
```

## Arguments

- se:

  A SummarizedExperiment object

- n:

  Number of top variable features to select (default: 500)

- assay_name:

  Name of assay to use (default: "counts")

## Value

A SummarizedExperiment subset to the top n variable features

## Examples

``` r
# Assuming 'se' is a SummarizedExperiment
data(example_se)
se_top <- top_variable_features(example_se, n = 500)
```
