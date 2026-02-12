
# Create example objects for each `from` class ----
gene_ids <- paste0("gene", seq_len(200))
cell_ids <- paste0("cell", seq_len(100))
mat <- matrix(rpois(20000, 5), ncol = 100, dimnames = list(gene_ids, cell_ids))

se <- SummarizedExperiment(assays = list(mat))
rse <- as(se, "RangedSummarizedExperiment")
sce <- as(se, "SingleCellExperiment")

# Start tests ----
test_that("coercion methods work", {
    
    ex1 <- as(se, "GraphExperiment")
    ex2 <- as(rse, "GraphExperiment")
    ex3 <- as(sce, "GraphExperiment")
    
    expect_s4_class(ex1, "GraphExperiment")
    expect_s4_class(ex2, "GraphExperiment")
    expect_s4_class(ex3, "GraphExperiment")
})
