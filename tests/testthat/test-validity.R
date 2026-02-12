
# Create example data ----
gene_ids <- paste0("gene", seq_len(200))
cell_ids <- paste0("cell", seq_len(100))
mat <- matrix(rpois(20000, 5), ncol = 100, dimnames = list(gene_ids, cell_ids))
g <- graph_from_adjacency_matrix(cor(t(mat)), weighted = TRUE)

# Start tests ----
test_that(".node_validity works", {
    
    mat2 <- mat[1:90, ]
    expect_error(
        GraphExperiment(assays = list(counts = mat2), graphs = list(cor = g))
    )
})

test_that(".graphs_validity works", {
    
    ge <- GraphExperiment(assays = list(counts = mat), graphs = list(cor = g))
    expect_error(
        graphs(ge) <- NA
    )
    expect_error(
        GraphExperiment(assays = list(counts = mat), graphs = list(cor = "graph"))
    )
    
})
