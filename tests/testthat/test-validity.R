
# Create example data ----
gene_ids <- paste0("gene", seq_len(200))
cell_ids <- paste0("cell", seq_len(100))
mat <- matrix(rpois(20000, 5), ncol = 100, dimnames = list(gene_ids, cell_ids))
g <- graph_from_adjacency_matrix(cor(t(mat)), weighted = TRUE)

# Start tests ----
test_that(".node_validity works", {
    
    mat2 <- mat[1:90, ]
    ## Nodes in graph but not in rownames
    expect_error(
        GraphExperiment(assays = list(counts = mat2), graphs = list(cor = g))
    )
    
    ## Nodes in rownames but not in graph
    mat3 <- mat[1:10, ]
    rownames(mat3) <- paste0("gene", LETTERS[1:10])
    mat3 <- rbind(mat, mat3)
    expect_error(
        GraphExperiment(assays = list(counts = mat3), graphs = list(cor = g))
    )
})

test_that(".graphs_validity works", {
    
    ge <- GraphExperiment(assays = list(counts = mat), graphs = list(cor = g))
    
    ## graphs slot is not list, SimpleList, or NULL
    expect_error(graphs(ge) <- NA)
    
    expect_error(
        GraphExperiment(assays = list(counts = mat), graphs = list(cor = "graph"))
    )
    
})
