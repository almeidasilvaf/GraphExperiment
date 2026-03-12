
#' @name GraphExperiment-subset
#' 
#' @title Subsetting `GraphExperiment` objects
#' 
#' @aliases [,GraphExperiment,ANY,ANY,ANY-method
#' 
#' @description
#' The subsetting method for \code{\link{GraphExperiment}} objects ensures
#' that nodes from `igraph` objects are filtered to match rows and columns 
#' with the remainder of the object.
#' 
#' @section subset:
#' \describe{
#' \item{\code{[}:}{ subsetting method}
#' }
#' 
#' @param x A \code{\link{GraphExperiment}} object.
#' @param i Numeric, row indices for subsetting.
#' @param j Numeric, column indices for subsetting.
#' 
#' @return a \code{\link{GraphExperiment}} object.
#' 
#' @examples
#' # Simulate elements of a GraphExperiment object
#' ## Assays
#' gene_ids <- paste0("gene", seq_len(200))
#' cell_ids <- paste0("cell", seq_len(100))
#' mat <- matrix(rpois(20000, 5), ncol = 100, dimnames = list(gene_ids, cell_ids))
#' 
#' ## Graph (with node attributes)
#' g <- graph_from_adjacency_matrix(cor(t(mat)), weighted = TRUE)
#' V(g)$degree <- igraph::strength(g)
#' 
#' ## rowData
#' rdata <- data.frame(
#'     row.names = gene_ids, 
#'     pathway = sample(c("P1", "P2"), size = length(gene_ids), replace = TRUE),
#'     coding = sample(c(TRUE, FALSE), size = length(gene_ids), replace = TRUE)
#' )
#' 
#' 
#' # Create a GraphExperiment object
#' ge <- GraphExperiment(
#'     assays = list(counts = mat), 
#'     rowData = rdata,
#'     graphs = list(cor = g)
#' )
#' ge
#' 
#' # Subset object
#' ge[1:5, 1:5]
NULL

#' @importFrom methods callNextMethod
#' @importFrom igraph induced_subgraph
#' @export
setMethod(
    "[", c("GraphExperiment", "ANY", "ANY"),
    function(x, i, j, ..., drop = FALSE) {
        
        if(missing(i)) i <- TRUE
        if(missing(j)) j <- TRUE
        
        glist <- graphs(x)
        x <- callNextMethod(x, i, j, ..., drop = drop)
        
        if(length(glist) > 0 || !isTRUE(i)) {
            keep <- rownames(x)
            fglist <- SimpleList(lapply(glist, function(g) {
                induced_subgraph(g, vids = keep)
            }))
            graphs(x) <- fglist
        }
        
        return(x)
    }
)
