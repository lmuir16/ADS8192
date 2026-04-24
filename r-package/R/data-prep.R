# R/data-prep.R

#' Select top variable features
#'
#' @param se A SummarizedExperiment object
#' @param n Number of top variable features to select (default: 500)
#' @param assay_name Name of assay to use (default: "counts")
#'
#' @return A SummarizedExperiment subset to the top n variable features
#'
#' @importFrom SummarizedExperiment assay assayNames
#' @importFrom stats var
#' @export
#'
#' @examples
#' # Assuming 'se' is a SummarizedExperiment
#' data(example_se)
#' se_top <- top_variable_features(example_se, n = 500)

top_variable_features <- function(se, n = 500, assay_name = "counts") {

  # Input checks
  if (!inherits(se, "SummarizedExperiment")) {
    stop("se must be a SummarizedExperiment object")
  }
  if (!assay_name %in% assayNames(se)) {
    stop(paste0("assay '", assay_name, "' not found in se"))
  }

  # Cap n at the number of available features
  if (n > nrow(se)) {
    warning(paste0("n (", n, ") exceeds the number of available features (",
                   nrow(se), "). Using all ", nrow(se), " features."))
    n <- nrow(se)
  }

  # Subset to top n variable genes (for tractable pairwise correlations)
  vars <- apply(assay(se, assay_name), 1, var)
  keep_genes <- names(sort(vars, decreasing = TRUE))[1:n]
  return(se[keep_genes, ])
}

#' Select set of random samples
#'
#' @param se A SummarizedExperiment object
#' @param n Number of random samples to select (default: 200)
#'
#' @return A SummarizedExperiment subset to n random samples
#'
#' @export
#'
#' @note For reproducibility, set a random seed with \code{set.seed()}
#'   before calling this function.

select_random_samples <- function(se, n = 200) {

  # Input checks
  if (!inherits(se, "SummarizedExperiment")) {
    stop("se must be a SummarizedExperiment object")
  }
  if (!is.numeric(n) || length(n) != 1 || n < 1) {
    stop("n must be a single positive integer")
  }

  # Subset to 200 random samples
  keep_samples <- sample(ncol(se), min(n, ncol(se)))
  se <- se[, keep_samples]
  return(se)
}


#' Prepare log-normalized expression matrix from a SummarizedExperiment
#'
#' Extracts the counts assay from a SummarizedExperiment object, applies
#' log2(x + 1) normalization, and remaps Ensembl row identifiers to HGNC
#' gene symbols using rowData. Falls back to Ensembl IDs for genes with
#' missing or empty symbols, and resolves duplicate symbols with
#' \code{make.unique()}.
#'
#' @param se A SummarizedExperiment object with a \code{counts} assay
#'   If \code{rowData()} contains a \code{gene_name} column, Ensembl IDs
#'   are remapped to HGNC gene symbols. Otherwise rownames are used as-is.
#'
#' @return A log2-normalized genes x samples matrix with gene symbols as
#'   row names
#'
#' @importFrom SummarizedExperiment assay assayNames rowData
#' @export
#'
#' @examples
#' library(SummarizedExperiment)
#' data(example_se)
#' mat <- prepare_expression_matrix(example_se)

prepare_expression_matrix <- function(se) {
  if (!inherits(se, "SummarizedExperiment")) {
    stop("se must be a SummarizedExperiment object")
  }
  if (!"counts" %in% assayNames(se)) {
    stop("se must contain a 'counts' assay")
  }
  # Log2-normalize counts
  mat <- log2(assay(se, "counts") + 1)
  # Remap row identifiers to gene symbols if gene_name is available
  # Falls back to existing rownames (e.g. Ensembl IDs) if not present
  if ("gene_name" %in% colnames(rowData(se))) {
    gene_symbols <- rowData(se)$gene_name
    gene_symbols[is.na(gene_symbols) | gene_symbols == ""] <-
      rownames(se)[is.na(gene_symbols) | gene_symbols == ""]
    rownames(mat) <- make.unique(gene_symbols)
  }
  return(mat)
}
