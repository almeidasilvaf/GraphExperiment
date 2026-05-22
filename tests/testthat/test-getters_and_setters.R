
# Construct `GraphExperiment` object ----
gene_ids <- paste0("gene", seq_len(200))
cell_ids <- paste0("cell", seq_len(100))
mat <- matrix(rpois(20000, 5), ncol = 100, dimnames = list(gene_ids, cell_ids))

g <- graph_from_adjacency_matrix(cor(t(mat)), weighted = TRUE)
V(g)$degree <- igraph::strength(g)

g2 <- graph_from_adjacency_matrix(cor(mat), weighted = TRUE)

rdata <- data.frame(
    row.names = gene_ids,
    pathway = sample(c("P1", "P2"), size = length(gene_ids), replace = TRUE),
    coding = sample(c(TRUE, FALSE), size = length(gene_ids), replace = TRUE)
)

cdata <- data.frame(
    row.names = cell_ids,
    celltype = sample(c("ct1", "ct2"), size = length(cell_ids), replace = TRUE)
)

ge <- GraphExperiment(
    assays = list(counts = mat),
    rowData = rdata,
    colData = cdata,
    rowGraphs = list(cor = g),
    colGraphs = list(cellcor = g2)
)
ge


# Start tests ----
test_that("Constructor function works", {
    ge1 <- GraphExperiment(
        assays = list(counts = mat), 
        rowGraphs = list(cor = g), colGraphs = list(cellcor = g2)
    )
    
    rg2 <- igraph::delete_vertices(g, paste0("gene", 1:10))
    cg2 <- igraph::delete_vertices(g2, paste0("cell", 1:10))
    
    expect_s4_class(ge1, "GraphExperiment")
    expect_error(
        GraphExperiment(assays = list(counts = mat), rowGraphs = "graph")
    )
    expect_error(
        GraphExperiment(assays = list(counts = mat), colGraphs = "graph")
    )
    
    expect_error(
        GraphExperiment(assays = list(counts = mat), rowGraphs = list("graph"))
    )
    expect_error(
        GraphExperiment(assays = list(counts = mat), colGraphs = list("graph"))
    )
    
    expect_error(
        GraphExperiment(assays = list(counts = mat), rowGraphs = rg2)
    )
    expect_error(
        GraphExperiment(assays = list(counts = mat), colGraphs = cg2)
    )
    
})

test_that("Getters work", {
    
    ## Basic GE object
    expect_true(is(rowGraphs(ge), "SimpleList"))
    expect_true(is(rowGraph(ge, 1), "igraph"))
    expect_true(is(rowGraph(ge), "igraph"))
    expect_equal(rowGraphNames(ge), "cor")
    
    expect_true(is(colGraphs(ge), "SimpleList"))
    expect_true(is(colGraph(ge, 1), "igraph"))
    expect_true(is(colGraph(ge), "igraph"))
    expect_equal(colGraphNames(ge), "cellcor")
    
    ## More complex GE objects
    ge2 <- ge
    V(rowGraph(ge2, 1))$rnum <- rnorm(vcount(rowGraph(ge2)))
    V(colGraph(ge2, 1))$cnum <- rnorm(vcount(colGraph(ge2)))
    
    expect_equal(length(rowData(ge)), 3)
    expect_true(length(rowData(ge2)) >3)
    
    at_names <- igraph::vertex_attr_names(rowGraph(ge))
    expect_true(all(gsub(".*__", "", colnames(rowData(ge))) %in% at_names))
    
    at_names2 <- igraph::vertex_attr_names(colGraph(ge))
    expect_true(all(gsub(".*__", "", colnames(colData(ge))) %in% at_names2))
    
    ## Empty graph
    ge_empty <- ge
    rowGraphs(ge_empty) <- NULL
    colGraphs(ge_empty) <- NULL
    
    expect_error(rowGraph(ge_empty))
    expect_error(colGraph(ge, 10))
})

test_that("Setters work", {
    
    nge <- ge
    rowGraphs(nge) <- list(g1 = g, g2 = g)
    rowGraphNames(nge) <- c("graph1", "graph2")
    rowGraph(nge, 1) <- g
    
    colGraphs(nge) <- list(g1 = g2, g2 = g2)
    colGraphNames(nge) <- c("cgraph1", "cgraph2")
    colGraph(nge, 1) <- g2
    
    expect_true(is(rowGraph(nge, 1), "igraph"))
    expect_true(is(rowGraphs(nge)[[1]], "igraph"))
    expect_equal(rowGraphNames(nge), c("graph1", "graph2"))
    expect_error(rowGraphs(nge) <- "graphs")
    expect_error(rowGraph(nge) <- "graph")
    expect_error(rowGraphNames(nge) <- c(1, 2))
    
    expect_true(is(colGraph(nge, 1), "igraph"))
    expect_true(is(colGraphs(nge)[[1]], "igraph"))
    expect_equal(colGraphNames(nge), c("cgraph1", "cgraph2"))
    expect_error(colGraphs(nge) <- "graphs")
    expect_error(colGraph(nge) <- "graph")
    expect_error(colGraphNames(nge) <- c(1, 2))
})

test_that("show method works", {

    expect_null(show(ge))
})

