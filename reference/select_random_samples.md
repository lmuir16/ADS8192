# Select set of random samples

Select set of random samples

## Usage

``` r
select_random_samples(se, n = 200)
```

## Arguments

- se:

  A SummarizedExperiment object

- n:

  Number of random samples to select (default: 200)

## Value

A SummarizedExperiment subset to n random samples

## Note

For reproducibility, set a random seed with
[`set.seed()`](https://rdrr.io/r/base/Random.html) before calling this
function.
