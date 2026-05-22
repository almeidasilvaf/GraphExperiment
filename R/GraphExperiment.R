
#' The `GraphExperiment` S4 class
#' 
#' The `GraphExperiment` class was designed to represent rectangular, 
#' quantitative data (e.g., from transcriptomics, proteomics, metabolomics) 
#' along with graphs showing how features (e.g., genes, proteins, compounds)
#' and observations (e.g., samples, cells, species, spots) interact with each 
#' other. It extends \code{SingleCellExperiment} by providing users with 
#' additional slots where row/column graphs can be stored.
#'
#' @param ... Arguments passed to the \code{SingleCellExperiment} constructor
#' function.
#' @param rowGraphs A list of `igraph` objects with one or multiple graphs
#' representing how features relate to each other. Node names (i.e., 
#' \code{V(graph)$name}) must match rownames of assays.
#' @param colGraphs A list of `igraph` objects with one or multiple graphs
#' representing how observations relate to each other. Node names (i.e., 
#' \code{V(graph)$name}) must match colnames of assays.
#' 
#' @return A \code{GraphExperiment} object.
#'
#' @details
#' Like \code{SingleCellExperiment}, the \code{GraphExperiment} S4 
#' class stores quantitative data with associated metadata (i.e., 
#' \code{rowData} and \code{colData}) along with embeddings from dimensionality 
#' reduction techniques. However, it provides users with additional 
#' containers for \code{igraph} objects containing graphs describing how 
#' features and/or observations interact with each other. Graphs for features 
#' are stored in a \code{rowGraphs} slot, and graphs for observations are
#' stored in a \code{colGraphs} slot. Both slots can hold one or 
#' multiple \code{igraph} objects with some sort of network representation of 
#' the features in \code{rownames} or observations in \code{colnames}. Example
#' graphs for features can be coexpression networks, regulatory networks, 
#' or co-abundance networks. Example graphs for observations can be 
#' cell-cell distances, sample-sample distances, or species networks 
#' representing genealogies (as an alternative to phylogenies).
#'
#' Importantly, node names in each row graph must match \code{rownames} of the
#' assays, and node names in each column graph must match \code{colnames} of
#' the assays. Besides, subsetting methods simultaneously subset \code{assays}, 
#' \code{rowData}, \code{rowGraphs}, \code{colData}, and \code{colGraphs}.
#' 
#' Besides the constructor function (\code{GraphExperiment()}), a 
#' \code{GraphExperiment} object can also be created by coercing from a
#' \code{SummarizedExperiment} or 
#' \code{SingleCellExperiment} object.
#' 
#' @rdname GraphExperiment
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom S4Vectors SimpleList
#' @importFrom methods new as
#' @export
#' @examples
#' # Example 1: from constructor function ----
#' ## Simulate a matrix with 200 genes and 100 cells
#' gene_ids <- paste0("gene", seq_len(200))
#' cell_ids <- paste0("cell", seq_len(100))
#' mat <- matrix(rpois(20000, 5), ncol = 100, dimnames = list(gene_ids, cell_ids))
#'
#' ## Create a rowGraph from correlations (`igraph` object)
#' g <- graph_from_adjacency_matrix(cor(t(mat)), weighted = TRUE)
#' 
#' ## Create a colGraph, also from correlations (but see scran::buildSNNGraph)
#' g2 <- graph_from_adjacency_matrix(cor(mat), weighted = TRUE)
#'
#' ## Construct `GraphExperiment` object
#' ge <- GraphExperiment(
#'     assays = list(counts = mat),
#'     rowGraphs = list(cor = g),
#'     colGraphs = list(cellcor = g2)
#' )
#' ge
#' 
#' # Example 2: From `SingleCellExperiment` object ----
#' sce <- SingleCellExperiment(assays = list(counts = mat))
#' ge <- as(sce, "GraphExperiment")
#' ge
#' 
GraphExperiment <- function(..., rowGraphs = list(), colGraphs = list()) {
    
    sce <- SingleCellExperiment(...)
    
    # rowGraphs
    if(is(rowGraphs, "list")) { rowGraphs <- SimpleList(rowGraphs) }
    if(is.null(names(rowGraphs)) & length(rowGraphs) >0) {
        names(rowGraphs) <- paste0("graph", seq_along(rowGraphs))
    }
    
    # colGraphs
    if(is(colGraphs, "list")) { colGraphs <- SimpleList(colGraphs) }
    if(is.null(names(colGraphs)) & length(colGraphs) >0) {
        names(colGraphs) <- paste0("graph", seq_along(colGraphs))
    }
    
    ge <- new(
        "GraphExperiment", sce, rowGraphs = rowGraphs, colGraphs = colGraphs
    )
    
    return(ge)
}

