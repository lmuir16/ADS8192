#' GTEx Skeletal Muscle SummarizedExperiment for testing
#'
#' A subset of the GTEx skeletal muscle bulk RNA-seq dataset accessed via
#' recount3. Contains the top 2000 most variable genes across 200 randomly
#' sampled donors, suitable for demonstrating gene co-expression network
#' analysis.
#'
#' @format A SummarizedExperiment with:
#' \describe{
#'   \item{assays}{counts - transformed count matrix (2000 genes x 200 samples)}
#'   \item{colData}{gtex.age, gtex.sex}
#'   \item{rowData}{gene_id, gene_name (HGNC symbol, with Ensembl ID fallback)}
#' }
#'
#' @source GTEx v8 skeletal muscle tissue via recount3
#'   (\url{https://rna.recount.bio/})
#'
#' @examples
#' library(SummarizedExperiment)
#' data(example_se)
#' example_se
#' colData(example_se)
"example_se"
