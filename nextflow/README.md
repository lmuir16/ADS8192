# Nextflow Pipeline

**corinet** provides a Nextflow workflow for co-expression network analysis from bulk RNA-seq data. It uses the corinet R package to compute pairwise gene correlations, detect modules, and produce network visualizations from a counts matrix, much like the corinet CLI implementation.

## Structure

- `main.nf` - Main workflow entry point
- `modules/` - Reusable process modules
- `data/` - Test data

## Setup

Install prerequisites:

1. [Docker](https://docs.docker.com/get-docker/)
2. [Nextflow](https://www.nextflow.io/docs/latest/install.html)

Build the Docker image:

```bash
docker build -t net-lmuir:0.0.1 ./r-package
docker run --rm net-lmuir:0.0.1 corinet --help
```

## Usage

Run the pipeline with test data:

```bash
nextflow run nextflow/main.nf -profile docker --counts nextflow/data/test/counts.tsv 
```

## Parameters

- `--counts` - Path to counts file (TSV format)

## Output

Results are written to `results/`:

- `gene_correlations.tsv` — gene pairs above correlation threshold
- `network_summary.tsv` — global network statistics
- `hub_genes.tsv` — top hub genes ranked by degree
- `module_assignments.tsv` — module membership for all genes
- `node_stats.tsv` — degree, betweenness, and module per gene
- `module_summary.tsv` — per-module size and top hub gene
- `correlation_heatmap.png` — gene-gene correlation heatmap
- `network_plot.png` — force-directed network graph