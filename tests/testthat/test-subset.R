
# Create example `GraphExperiment` object
gene_ids <- paste0("gene", seq_len(200))
cell_ids <- paste0("cell", seq_len(100))
mat <- matrix(rpois(20000, 5), ncol = 100, dimnames = list(gene_ids, cell_ids))

## Graph (with node attributes)
g <- graph_from_adjacency_matrix(cor(t(mat)), weighted = TRUE)
V(g)$degree <- igraph::strength(g)

## rowData
rdata <- data.frame(
    row.names = gene_ids,
    pathway = sample(c("P1", "P2"), size = length(gene_ids), replace = TRUE),
    coding = sample(c(TRUE, FALSE), size = length(gene_ids), replace = TRUE)
)


# Create a GraphExperiment object
ge <- GraphExperiment(
    assays = list(counts = mat),
    rowData = rdata,
    graphs = list(cor = g)
)
ge

# Start tests ----
test_that("Subsetting rows also subsets nodes", {
    
    fge <- ge[1:5, 1:5]
    
    expect_equal(nrow(fge), 5)
    expect_equal(vcount(graph(fge, 1)), 5)
    
    expect_equal(dim(ge[, ]), dim(ge))
})
