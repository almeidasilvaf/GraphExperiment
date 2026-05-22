
#' Methods for `GraphExperiment` objects
#' 
#' The \code{GraphExperiment} class provides users with methods to get and
#' set graphs (\code{igraph} objects) representing how features and observations 
#' of \code{SingleCellExperiment} objects relate to each other.
#'
#' @param x A \code{GraphExperiment} object.
#' @param i List element (numeric for index, character for name) of the element
#' to access or replace.
#' @param value Replacement value for replacement methods.
#' @param ... Ignored.
#' @return Return values depend on the method. See details and examples.
#' 
#' @section rowGraphs/colGraphs and rowGraph/colGraph methods:
#' \describe{
#'   \item{\code{rowGraphs(x)} and \code{colGraphs(x)}: }{
#'     Getter for a \code{SimpleList} of \code{igraph} objects representing 
#'     rows and columns, respectively.}
#'   \item{\code{rowGraphs(x) <- value} and \code{colGraphs(x) <- value}: }{
#'     Setter for a SimpleList or list (coerced to SimpleList) of \code{igraph}
#'     objects representing rows and columns, respectively.}
#'   \item{\code{rowGraph(x, i)} and \code{colGraph(x, i)}: }{
#'     Getter for an \code{igraph} object containing graph `i` from the list
#'     stored in \code{rowGraphs} and \code{colGraphs}, respectively.}
#'   \item{\code{rowGraph(x, i) <- value} and \code{colGraph(x, i) <- value}: }{
#'     Setter for an \code{igraph} object to be stored in element `i` of 
#'     \code{rowGraphs} and \code{colGraphs}.}
#' }
#' 
#' @section rowGraphNames and colGraphNames methods:
#' \describe{
#'   \item{\code{rowGraphNames(x)} and \code{colGraphNames(x)}: }{
#'     Getter to extract names of graphs in \code{rowGraphs} 
#'     and \code{colGraphs}.}
#'   \item{\code{rowGraphNames(x) <- value} and \code{colGraphNames(x) <- value}: }{
#'     Setter to assign new names to the graphs stored in \code{rowGraphs}
#'     and \code{colGraphs}.}
#' }
#' 
#' @section rowData method:
#' \describe{
#'   \item{\code{rowData(x)}: }{
#'     Getter to extract rowData (as in \code{SingleCellExperiment} objects),
#'     but with node attributes of graphs included.}
#' }
#' 
#' @section colData method:
#' \describe{
#'   \item{\code{colData(x)}: }{
#'     Getter to extract colData (as in \code{SingleCellExperiment} objects),
#'     but with node attributes of graphs included.}
#' }
#' 
#' @name GraphExperiment-methods
#' @aliases 
#' rowGraphs rowGraphs<- colGraphs colGraphs<-
#' rowGraph rowGraph<- colGraph colGraph<-
#' rowGraphNames rowGraphNames<- colGraphNames colGraphNames<-
#' rowData colData
#' 
#' @examples 
#' # Simulate elements of a GraphExperiment object
#' ## Assays
#' gene_ids <- paste0("gene", seq_len(200))
#' cell_ids <- paste0("cell", seq_len(100))
#' mat <- matrix(rpois(20000, 5), ncol = 100, dimnames = list(gene_ids, cell_ids))
#' 
#' ## rowGraph (with node attributes)
#' g <- graph_from_adjacency_matrix(cor(t(mat)), weighted = TRUE)
#' V(g)$degree <- igraph::strength(g)
#' 
#' ## colGraph
#' g2 <- graph_from_adjacency_matrix(cor(mat), weighted = TRUE)
#' 
#' ## rowData
#' rdata <- data.frame(
#'     row.names = gene_ids, 
#'     pathway = sample(c("P1", "P2"), size = length(gene_ids), replace = TRUE),
#'     coding = sample(c(TRUE, FALSE), size = length(gene_ids), replace = TRUE)
#' )
#' 
#' ## colData
#' cdata <- data.frame(
#'     row.names = cell_ids, 
#'     celltype = sample(c("ct1", "ct2"), size = length(cell_ids), replace = TRUE)
#' )
#' 
#' 
#' # Create a GraphExperiment object
#' ge <- GraphExperiment(
#'     assays = list(counts = mat), 
#'     rowData = rdata,
#'     colData = cdata,
#'     rowGraphs = list(cor = g),
#'     colGraphs = list(cellcor = g2)
#' )
#' ge
#' 
#' # Extract graph names
#' rowGraphNames(ge)
#' colGraphNames(ge)
#' 
#' 
#' # Extract graphs
#' rowGraphs(ge)
#' rowGraph(ge, "cor")
#' 
#' colGraphs(ge)
#' colGraph(ge, "cellcor")
#' 
#' 
#' # Add a new graph
#' rowGraph(ge, "newcor") <- g
#' ge
#' 
#' # Add a list of graphs
#' colGraphs(ge) <- list(cellcor = g2, new_cellcor = g2)
#' ge
#' 
#' # Replace graph names
#' rowGraphNames(ge) <- c("network1", "network2")
#' ge
#' 
#' # Access rowData (note: rowData + node attributes combined)
#' rowData(ge)
#' 
#' # Access colData (note: colData + node attributes combined)
#' colData(ge)
#'
NULL

## Getters ---------------------------------------------------------------------


# rowData and colData ----------------------------------------------------------

#' @param use.names Passed to the \code{rowData} method of
#' \code{SingleCellExperiment}. Default: TRUE.
#' @rdname GraphExperiment-methods
#' @importFrom SummarizedExperiment rowData
#' @export
setMethod(
    "rowData", "GraphExperiment",
    function(x, use.names = TRUE, ...) {
        rdata <- callNextMethod()
        rdata <- .get_metadata(x@rowGraphs, rdata)
        
        return(rdata)
    }
)

#' @param use.names Passed to the \code{colData} method of
#' \code{SingleCellExperiment}. Default: TRUE.
#' @rdname GraphExperiment-methods
#' @importFrom SummarizedExperiment colData
#' @export
setMethod(
    "colData", "GraphExperiment",
    function(x, use.names = TRUE, ...) {
        cdata <- callNextMethod()
        cdata <- .get_metadata(x@colGraphs, cdata)
        
        return(cdata)
    }
)



# rowGraphs and colGraphs ------------------------------------------------------

#' @rdname GraphExperiment-methods
#' @export
setMethod(
    "rowGraphs", "GraphExperiment", 
    function(x) {
        glist <- x@rowGraphs
        if(length(glist) >0) {
            rdata <- rowData(x)
            glist <- SimpleList(lapply(glist, function(x) {
                .metadata2nat(rdata, x)
            }))
        }
        return(glist)
    }
)

#' @rdname GraphExperiment-methods
#' @export
setMethod(
    "colGraphs", "GraphExperiment", 
    function(x) {
        glist <- x@colGraphs
        if(length(glist) >0) {
            cdata <- colData(x)
            glist <- SimpleList(lapply(glist, function(x) {
                .metadata2nat(cdata, x)
            }))
        }
        return(glist)
    }
)


# rowGraph and colGraph --------------------------------------------------------

#' @rdname GraphExperiment-methods
#' @export
setMethod(
    "rowGraph", c("GraphExperiment", "missing"),
    function(x, i) {
        
        glist <- x@rowGraphs
        if(length(glist) == 0) {
            stop("The 'rowGraphs' slot is empty.")
        }
        .metadata2nat(rowData(x), glist[[1]])
    }
)

#' @rdname GraphExperiment-methods
#' @export
setMethod(
    "rowGraph", "GraphExperiment", 
    function(x, i) {
        
        tryCatch({
            .metadata2nat(rowData(x), x@rowGraphs[[i]])
        }, error = function(e) {
            stop("Invalid index or name. Could not find element ", i, " in 'rowGraphs'.")
        })
    } 
)

#' @rdname GraphExperiment-methods
#' @export
setMethod(
    "colGraph", c("GraphExperiment", "missing"),
    function(x, i) {
        
        glist <- x@colGraphs
        if(length(glist) == 0) {
            stop("The 'colGraphs' slot is empty.")
        }
        .metadata2nat(colData(x), glist[[1]])
    }
)

#' @rdname GraphExperiment-methods
#' @export
setMethod(
    "colGraph", "GraphExperiment", 
    function(x, i) {
        
        tryCatch({
            .metadata2nat(colData(x), x@colGraphs[[i]])
        }, error = function(e) {
            stop("Invalid index or name. Could not find element ", i, " in 'colGraphs'.")
        })
    } 
)

# rowGraphNames and colGraphNames ----------------------------------------------

#' @rdname GraphExperiment-methods
#' @export
setMethod(
    "rowGraphNames", "GraphExperiment",
    function(x) {
        names(rowGraphs(x))
    }
)

#' @rdname GraphExperiment-methods
#' @export
setMethod(
    "colGraphNames", "GraphExperiment",
    function(x) {
        names(colGraphs(x))
    }
)


## Setters ----

#' @rdname GraphExperiment-methods
#' @export
#' @importFrom BiocBaseUtils setSlots
setReplaceMethod(
    "rowGraphs", "GraphExperiment",
    function(x, value) {
        
        if(is(value, "list")) { value <- SimpleList(value) }
        setSlots(x, rowGraphs = value)
    }
)

#' @rdname GraphExperiment-methods
#' @export
#' @importFrom BiocBaseUtils setSlots
setReplaceMethod(
    "colGraphs", "GraphExperiment",
    function(x, value) {
        
        if(is(value, "list")) { value <- SimpleList(value) }
        setSlots(x, colGraphs = value)
    }
)




#' @rdname GraphExperiment-methods
#' @export
setReplaceMethod(
    "rowGraph", "GraphExperiment",
    function(x, i, value) {
        rowGraphs(x)[[i]] <- value
        x
    }
)

#' @rdname GraphExperiment-methods
#' @export
setReplaceMethod(
    "colGraph", "GraphExperiment",
    function(x, i, value) {
        colGraphs(x)[[i]] <- value
        x
    }
)




#' @rdname GraphExperiment-methods
#' @export
setReplaceMethod(
    "rowGraphNames", c("GraphExperiment", "character"),
    function(x, value) {
        names(rowGraphs(x)) <- value
        x
    }
)

#' @rdname GraphExperiment-methods
#' @export
setReplaceMethod(
    "colGraphNames", c("GraphExperiment", "character"),
    function(x, value) {
        names(colGraphs(x)) <- value
        x
    }
)

