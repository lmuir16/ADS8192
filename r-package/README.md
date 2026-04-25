
<!-- README.md is generated from README.Rmd. Please edit that file -->

# corinet

<!-- badges: start -->

[![R-CMD-check](https://github.com/lmuir16/corinet/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/lmuir16/corinet/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/lmuir16/corinet/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/lmuir16/corinet/actions/workflows/pkgdown.yaml)
<!-- badges: end -->

**corinet** provides a lightweight pipeline for building **cor**relation
**net**works to analyze gene co-expression in bulk RNA-seq data.
Starting from a `SummarizedExperiment` object, the package supports the
full workflow from data preparation through network construction, module
detection, visualization, and export.

## Installation

corinet requires several Bioconductor dependencies. Install them first:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("SummarizedExperiment", "recount3"))
```

Then install corinet from GitHub:

``` r
# Using pak (recommended)
install.packages("pak")
pak::pak("lmuir16/corinet")

# Or using devtools
install.packages("devtools")
devtools::install_github("lmuir16/corinet")
```

## Quick Start

``` r
library(corinet)

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

For non-interactive use, corinet provides a CLI via Rapp. After
installing the package, install the CLI with:

``` r
Rapp::install_pkg_cli_apps("corinet")
```

All subcommands accept a genes x samples counts matrix as TSV or CSV
input:

``` bash
# Full network analysis pipeline — exports 6 TSV files
corinet network --counts counts.tsv --output results/ --n-top 500

# Correlation heatmap
corinet heatmap --counts counts.tsv --output results/ --n-top 500

# Force-directed network graph
corinet plot-network --counts counts.tsv --output results/ --n-top 500 --n-plot 80
```

See `corinet --help` and `corinet <subcommand> --help` for all available
options.

A small example counts matrix is available at
`tests/testdata/counts.tsv` in the package repository for testing
purposes.

## Learn More

See `vignette("getting-started", package = "corinet")` for a full
walkthrough of the pipeline.
