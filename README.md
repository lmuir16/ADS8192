
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ADS8192

<!-- badges: start -->

[![R-CMD-check](https://github.com/lmuir16/ADS8192/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/lmuir16/ADS8192/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/lmuir16/ADS8192/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/lmuir16/ADS8192/actions/workflows/pkgdown.yaml)
<!-- badges: end -->

ADS8192 provides a lightweight pipeline for building and analyzing gene
co-expression networks from bulk RNA-seq data. Starting from a
`SummarizedExperiment` object, the package supports the full workflow
from data preparation through network construction, module detection,
visualization, and export.

## Installation

ADS8192 requires several Bioconductor dependencies. Install them first:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("SummarizedExperiment", "recount3"))
```

Then install ADS8192 from GitHub:

``` r
# Using pak (recommended)
install.packages("pak")
pak::pak("lmuir16/ADS8192")

# Or using devtools
install.packages("devtools")
devtools::install_github("lmuir16/ADS8192")
```

## Quick Start

``` r
library(ADS8192)

# Load example GTEx skeletal muscle data
data(example_se)

# Use a small subset for a quick demonstration
se_top <- top_variable_features(example_se[1:100, 1:50], n = 50)
mat    <- prepare_expression_matrix(se_top)

# Build network
cor_mat <- pairwise_corr(mat)
adj_mat <- build_adjacency_mat(cor_mat, cor_threshold = 0.3)
net     <- build_network(adj_mat)

# Detect modules and summarize
gene_mods   <- detect_network_modules(net)
#> Modularity: 0.211
#> Modules found: 3
net_summary <- summarize_network(net, gene_mods$modules, gene_mods$modularity)
#>          metric    value
#> 1       n_genes  50.0000
#> 2       n_edges 598.0000
#> 3   mean_degree  23.9200
#> 4 median_degree  26.0000
#> 5  n_hubs_deg20  34.0000
#> 6     n_modules   3.0000
#> 7    modularity   0.2107

# Identify hub genes
hub_genes <- get_hub_genes(net, modules = gene_mods$modules, n = 10)
head(hub_genes)
#>        gene degree betweenness module
#> TPM2   TPM2     40    46.30589     M3
#> CKM     CKM     36    30.41338     M3
#> MB       MB     34    20.32788     M2
#> FHL1   FHL1     34    41.02601     M3
#> TTN     TTN     33    20.73674     M1
#> ACTA1 ACTA1     33    15.92916     M2
```

## Command Line Interface

For non-interactive use, ADS8192 provides a CLI via Rapp. After
installing the package, install the CLI with:

``` r
Rapp::install_pkg_cli_apps("ADS8192")
```

Example test data is included at `tests/testdata/counts.tsv`. Run the
pipeline from the terminal:

All subcommands accept a genes x samples counts matrix as input:

``` bash
# Full network analysis pipeline — exports 6 TSV files
ADS8192 network --counts tests/testdata/counts.tsv --output results/ --n-top 50

# Correlation heatmap
ADS8192 heatmap --counts tests/testdata/counts.tsv --output results/ --n-top 50

# Force-directed network graph
ADS8192 plot-network --counts tests/testdata/counts.tsv --output results/ --n-top 50 --n-plot 20
```

See `ADS8192 --help` and `ADS8192 <subcommand> --help` for all available
options.

## Learn More

See `vignette("getting-started", package = "ADS8192")` for a full
walkthrough of the pipeline.
