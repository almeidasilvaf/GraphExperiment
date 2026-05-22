
# Create example data ----
gene_ids <- paste0("gene", seq_len(200))
cell_ids <- paste0("cell", seq_len(100))
mat <- matrix(rpois(20000, 5), ncol = 100, dimnames = list(gene_ids, cell_ids))
rg <- graph_from_adjacency_matrix(cor(t(mat)), weighted = TRUE)
cg <- graph_from_adjacency_matrix(cor(mat), weighted = TRUE)

# Start tests ----
test_that(".node_validity works", {
    
    ## Nodes in graph but not in rownames
    expect_error(
        GraphExperiment(assays = list(counts = mat[1:90, ]), rowGraphs = list(cor = rg))
    )
    
    ## Nodes in graph but not in colnames
    expect_error(
        GraphExperiment(assays = list(counts = mat[, 1:90]), colGraphs = list(cor = cg))
    )
    
    ## Nodes in rownames but not in graph
    mat3 <- mat[1:10, ]
    rownames(mat3) <- paste0("gene", LETTERS[1:10])
    mat3 <- rbind(mat, mat3)
    expect_error(
        GraphExperiment(assays = list(counts = mat3), rowGraphs = list(cor = rg))
    )
    
    ## Nodes in colnames but not in graph
    mat4 <- mat[, 1:10]
    colnames(mat4) <- paste0("cell", LETTERS[1:10])
    mat4 <- cbind(mat, mat4)
    expect_error(
        GraphExperiment(assays = list(counts = mat4), colGraphs = list(cor = cg))
    )
    
})

test_that(".graphs_validity works", {
    
    ge <- GraphExperiment(assays = list(counts = mat), rowGraphs = list(cor = rg))
    
    ## graphs slot is not list, SimpleList, or NULL
    expect_error(rowGraphs(ge) <- NA)
    expect_error(colGraphs(ge) <- NA)
    
    expect_error(
        GraphExperiment(assays = list(counts = mat), rowGraphs = list(cor = "graph"))
    )
    
})
