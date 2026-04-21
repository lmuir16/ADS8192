# ADS8192

ADS8192 provides a lightweight pipeline for building and analyzing gene
co-expression networks from bulk RNA-seq data. Starting from a
`SummarizedExperiment` object, the package supports the full workflow
from data preparation through network construction, module detection,
visualization, and export.

## Features

- Select top variable genes for network construction
- Compute pairwise Pearson or Spearman gene-gene correlations
- Threshold correlations into a binary adjacency matrix
- Detect gene modules using Louvain or walktrap community detection
- Visualize networks as correlation heatmaps or force-directed graphs
- Export results to TSV files for downstream analysis

## Installation

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("SummarizedExperiment", "recount3"))

pak::pak("lmuir16/ADS8192")
```

## Getting Started

See the [Getting Started vignette](articles/getting-started.html) for a 
full walkthrough of the pipeline using GTEx skeletal muscle data.
