
#' Methods for `GraphExperiment` objects
#' 
#' The \code{GraphExperiment} class provides users with methods to get and
#' set graphs (\code{igraph} objects) representing how features of 
#' \linkS4class{SingleCellExperiment} objects relate to each other.
#'
#' @param x A \code{GraphExperiment} object.
#' @param i List element (numeric for index, character for name) of the element
#' to access or replace.
#' @param value Replacement value for replacement methods
#' 
#' @return Return values depend on the method. See details and examples.
#' 
#' @section graphs and graph methods:
#' \describe{
#'   \item{\code{graphs(x)}: }{
#'     Getter for a \code{SimpleList} of \code{igraph} objects.}
#'   \item{\code{graphs(x) <- value}: }{
#'     Setter for a SimpleList or list (coerced to SimpleList) of \code{igraph}
#'     objects.}
#'   \item{\code{graph(x, i)}: }{
#'     Getter for an \code{igraph} object containing graph `i` from the list
#'     stored in \code{graphs}.}
#'   \item{\code{graph(x, i) <- value}: }{
#'     Setter for an \code{igraph} object to be stored in element `i` of the 
#'     list \code{graphs}}
#' }
#' 
#' @section graphNames methods:
#' \describe{
#'   \item{\code{graphNames(x)}: }{
#'     Getter to extract names of graphs in the \code{graphs} slot.}
#'   \item{\code{graphNames(x) <-  value}: }{
#'     Setter to assign new names to the graphs stored in the \code{graphs}
#'     slot.}
#' }
#' 
#' @name GraphExperiment-methods
#' @aliases 
#' graphs graphs<-
#' graph graph<-
#' graphNames graphNames<-
#' 
#' 
#' @examples 
#' # Create a GraphExperiment object
#' gene_ids <- paste0("gene", seq_len(200))
#' cell_ids <- paste0("cell", seq_len(100))
#' mat <- matrix(rpois(20000, 5), ncol = 100, dimnames = list(gene_ids, cell_ids))
#' g <- graph_from_adjacency_matrix(cor(t(mat)), weighted = TRUE)
#' 
#' ge <- GraphExperiment(assays = list(counts = mat), graphs = list(cor = g))
#' ge
#' 
#' # Extract graph names
#' graphNames(ge)
#' ge
#' 
#' # Extract graphs
#' graphs(ge)
#' graph(ge, "cor")
#' 
#' # Add a new graph
#' graph(ge, "newcor") <- g
#' ge
#' 
#' # Add a list of graphs
#' graphs(ge) <- list(cor = g, newcor = g)
#' ge
#' 
#' # Replace graph names
#' graphNames(ge) <- c("network1", "network2")
#' ge
#'
NULL

## Getters ---------------------------------------------------------------------

#' @rdname GraphExperiment-methods
#' @export
setMethod(
    "graphs", "GraphExperiment", 
    function(x) x@graphs
)

#' @rdname GraphExperiment-methods
#' @export
setMethod(
    "graph", "GraphExperiment", 
    function(x, i) {
        
        tryCatch({
            graphs(x)[[i]]
        }, error = function(e) {
            stop("Invalid index or name. Could not find element ", i, " in list of graphs.")
        })
    } 
)


#' @rdname GraphExperiment-methods
#' @export
setMethod(
    "graphNames", "GraphExperiment",
    function(x) {
        names(graphs(x))
    }
)

## Setters ----

#' @rdname GraphExperiment-methods
#' @export
#' @importFrom BiocBaseUtils setSlots
setReplaceMethod(
    "graphs", "GraphExperiment",
    function(x, value) {
        
        if(is(value, "list")) { value <- SimpleList(value) }
        setSlots(x, graphs = value)
    }
)

#' @rdname GraphExperiment-methods
#' @export
setReplaceMethod(
    "graph", "GraphExperiment",
    function(x, i, value) {
        graphs(x)[[i]] <- value
        x
    }
)

#' @rdname GraphExperiment-methods
#' @export
setReplaceMethod(
    "graphNames", c("GraphExperiment", "character"),
    function(x, value) {
        names(graphs(x)) <- value
        x
    }
)
