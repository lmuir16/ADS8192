# GTEx Skeletal Muscle SummarizedExperiment for testing

A subset of the GTEx skeletal muscle bulk RNA-seq dataset accessed via
recount3. Contains the top 2000 most variable genes across 200 randomly
sampled donors, suitable for demonstrating gene co-expression network
analysis.

## Usage

``` r
example_se
```

## Format

A SummarizedExperiment with:

- assays:

  counts - transformed count matrix (2000 genes x 200 samples)

- colData:

  gtex.age, gtex.sex

- rowData:

  gene_id, gene_name (HGNC symbol, with Ensembl ID fallback)

## Source

GTEx v8 skeletal muscle tissue via recount3 (<https://rna.recount.bio/>)

## Examples

``` r
library(SummarizedExperiment)
#> Loading required package: MatrixGenerics
#> Loading required package: matrixStats
#> 
#> Attaching package: ‘MatrixGenerics’
#> The following objects are masked from ‘package:matrixStats’:
#> 
#>     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
#>     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
#>     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
#>     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
#>     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
#>     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
#>     colWeightedMeans, colWeightedMedians, colWeightedSds,
#>     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
#>     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
#>     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
#>     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
#>     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
#>     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
#>     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
#>     rowWeightedSds, rowWeightedVars
#> Loading required package: GenomicRanges
#> Loading required package: stats4
#> Loading required package: BiocGenerics
#> Loading required package: generics
#> 
#> Attaching package: ‘generics’
#> The following objects are masked from ‘package:base’:
#> 
#>     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
#>     setequal, union
#> 
#> Attaching package: ‘BiocGenerics’
#> The following objects are masked from ‘package:stats’:
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from ‘package:base’:
#> 
#>     Filter, Find, Map, Position, Reduce, anyDuplicated, aperm, append,
#>     as.data.frame, basename, cbind, colnames, dirname, do.call,
#>     duplicated, eval, evalq, get, grep, grepl, is.unsorted, lapply,
#>     mapply, match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
#>     rank, rbind, rownames, sapply, saveRDS, table, tapply, unique,
#>     unsplit, which.max, which.min
#> Loading required package: S4Vectors
#> 
#> Attaching package: ‘S4Vectors’
#> The following object is masked from ‘package:utils’:
#> 
#>     findMatches
#> The following objects are masked from ‘package:base’:
#> 
#>     I, expand.grid, unname
#> Loading required package: IRanges
#> Loading required package: Seqinfo
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> 
#> Attaching package: ‘Biobase’
#> The following object is masked from ‘package:MatrixGenerics’:
#> 
#>     rowMedians
#> The following objects are masked from ‘package:matrixStats’:
#> 
#>     anyMissing, rowMedians
data(example_se)
example_se
#> class: RangedSummarizedExperiment 
#> dim: 2000 200 
#> metadata(8): time_created recount3_version ... annotation recount3_url
#> assays(1): counts
#> rownames(2000): ENSG00000143632.14 ENSG00000198804.2 ...
#>   ENSG00000128578.9 ENSG00000174903.15
#> rowData names(10): source type ... havana_gene tag
#> colnames(200): GTEX-S3XE-2026-SM-3K2B5.1 GTEX-1MJIX-2126-SM-E9J3A.1 ...
#>   GTEX-1E2YA-0526-SM-7P8RY.1 GTEX-S3LF-0526-SM-EZ6LK.1
#> colData names(2): gtex.age gtex.sex
colData(example_se)
#> DataFrame with 200 rows and 2 columns
#>                               gtex.age  gtex.sex
#>                            <character> <integer>
#> GTEX-S3XE-2026-SM-3K2B5.1        50-59         1
#> GTEX-1MJIX-2126-SM-E9J3A.1       40-49         1
#> GTEX-QV44-2026-SM-2S1RD.1        50-59         1
#> GTEX-14ASI-0526-SM-5QGQP.1       60-69         1
#> GTEX-SNMC-1426-SM-2XCFM.1        20-29         1
#> ...                                ...       ...
#> GTEX-1HSKV-0126-SM-E6CI1.1       60-69         1
#> GTEX-13CF3-1826-SM-5J1NK.1       60-69         2
#> GTEX-16GPK-0526-SM-7KUER.1       60-69         1
#> GTEX-1E2YA-0526-SM-7P8RY.1       50-59         1
#> GTEX-S3LF-0526-SM-EZ6LK.1        70-79         1
```
