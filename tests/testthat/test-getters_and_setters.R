
# Construct `GraphExperiment` object ----
gene_ids <- paste0("gene", seq_len(200))
cell_ids <- paste0("cell", seq_len(100))
mat <- matrix(rpois(20000, 5), ncol = 100, dimnames = list(gene_ids, cell_ids))
g <- graph_from_adjacency_matrix(cor(t(mat)), weighted = TRUE)
 
ge <- GraphExperiment(assays = list(counts = mat), graphs = list(cor = g))

rdata <- data.frame(
    row.names = gene_ids, 
    pathway = sample(c("P1", "P2"), size = length(gene_ids), replace = TRUE),
    coding = sample(c(TRUE, FALSE), size = length(gene_ids), replace = TRUE)
)


# Start tests ----
test_that("Constructor function works", {
    ge1 <- GraphExperiment(assays = list(counts = mat), graphs = list(cor = g))
    g2 <- igraph::delete_vertices(g, paste0("gene", 1:10))
    
    expect_s4_class(ge1, "GraphExperiment")
    expect_error(
        GraphExperiment(assays = list(counts = mat), graphs = "graph")
    )
    expect_error(
        GraphExperiment(assays = list(counts = mat), graphs = list("graph"))
    )
    expect_error(
        GraphExperiment(assays = list(counts = mat), graphs = g2)
    )
})

test_that("Getters work", {
    
    ## Basic GE object
    expect_true(is(graphs(ge), "SimpleList"))
    expect_true(is(graph(ge, 1), "igraph"))
    expect_true(is(graph(ge), "igraph"))
    expect_equal(graphNames(ge), "cor")
    
    ## More complex GE objects
    ge2 <- ge
    V(graph(ge2, 1))$degree <- degree(graph(ge2, 1))
    
    ge3 <- GraphExperiment(
        assays = list(counts = mat), 
        rowData = rdata,
        graphs = list(cor = g)
    )
    
    expect_equal(length(rowData(ge)), 0)
    expect_equal(length(rowData(ge2)), 1)
    
    at_names <- igraph::vertex_attr_names(graph(ge3))
    expect_true(all(colnames(rowData(ge3)) %in% at_names))
    
    ## Empty graph
    ge_empty <- ge
    graphs(ge_empty) <- NULL
    
    expect_error(graph(ge_empty))
    expect_error(graph(ge, 10))
})

test_that("Setters work", {
    
    nge <- ge
    graphs(nge) <- list(g1 = g, g2 = g)
    graphNames(nge) <- c("graph1", "graph2")
    graph(nge, 1) <- g
    
    expect_true(is(graph(nge, 1), "igraph"))
    expect_true(is(graphs(nge)[[1]], "igraph"))
    expect_equal(graphNames(nge), c("graph1", "graph2"))
    expect_error(graphs(nge) <- "graphs")
    expect_error(graph(nge) <- "graph")
    expect_error(graphNames(nge) <- c(1, 2))
})

test_that("show method works", {

    expect_null(show(ge))
})

