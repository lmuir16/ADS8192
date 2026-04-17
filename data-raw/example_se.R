## data-raw/example_se.R

BiocManager::install("recount3")
library(recount3)
library(SummarizedExperiment)

# Access GTEx skeletal muscle data
human_projects <- available_projects()
proj <- subset(human_projects, file_source == "gtex" &
                 project == "MUSCLE")
rse <- create_rse(proj)

# Transform raw counts
assay(rse, "counts") <- transform_counts(rse)

# Subset to top 2000 genes by variance for tractable network computation
vars <- apply(assay(rse, "counts"), 1, var)
keep_genes <- names(sort(vars, decreasing = TRUE))[1:2000]

# Subset to 200 samples for package size
set.seed(42)
keep_samples <- sample(ncol(rse), min(200, ncol(rse)))

example_se <- rse[keep_genes, keep_samples]
assays(example_se) <- list(counts = assay(example_se, "counts"))
colData(example_se) <- colData(example_se)[, c("gtex.age", "gtex.sex"),
                                           drop = FALSE]

usethis::use_data(example_se, overwrite = TRUE)
