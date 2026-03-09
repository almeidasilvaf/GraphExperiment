
#' Methods for `GraphExperiment` objects
#' 
#' The \code{GraphExperiment} class provides users with methods to get and
#' set graphs (\code{igraph} objects) representing how features of 
#' \code{\linkS4class{SingleCellExperiment}} objects relate to 
#' each other.
#'
#' @param x A \code{GraphExperiment} object.
#' @param i List element (numeric for index, character for name) of the element
#' to access or replace.
#' @param value Replacement value for replacement methods
#' @param ... Ignored.
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
#' rowData
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
#' # Access rowData (note: rowData + node attributes combined)
#' rowData(ge)
#'
NULL

## Getters ---------------------------------------------------------------------

#' Add rowData variables to node attributes
#' @importFrom igraph set_vertex_attr vertex_attr_names
#' @noRd
.rowdata2nat <- function(rowdata, graph) {
    
    # Remove rowData vars extracted from graph node attributes
    idx_remove <- grepl("__", names(rowdata))
    if(sum(idx_remove) >0) { rowdata <- rowdata[, !idx_remove, drop = FALSE] }
    
    # Add missing rowData vars to node attributes
    fgraph <- graph
    if(length(rowdata) >0) {
        nat_list <- as.list(rowdata)
        
        if(length(nat_list) >0) {
            idx <- match(rownames(rowdata), V(graph)$name)
            for(rvar in names(nat_list)) {
                fgraph <- set_vertex_attr(
                    fgraph, name = rvar, index = idx, value = nat_list[[rvar]]
                )
            }
        }
    }
    
    return(fgraph)
}

#' @param use.names Passed to the \code{rowData} method of
#' \code{\linkS4class{SingleCellExperiment}}. Default: TRUE.
#' @rdname GraphExperiment-methods
#' @importFrom SummarizedExperiment rowData
#' @export
setMethod(
    "rowData", "GraphExperiment",
    function(x, use.names = TRUE, ...) {
        rdata <- callNextMethod()
        
        glist <- x@graphs
        gnames <- names(glist)
        
        # Get a data frame of node attributes for all graphs (if any)
        nodeat_df <- NULL
        if(length(glist) > 0) {
            nodeat_df <- Reduce(cbind, lapply(seq_along(glist), function(n) {
                df <- igraph::as_data_frame(glist[[n]], what = "vertices")
                names(df) <- paste0(gnames[n], "__", names(df))
                df[[1]] <- NULL
                return(df)
            }))
        }
        
        # Extend rowData to include node attributes
        if(length(nodeat_df) > 0) {
            idx_cbind <- match(rownames(rdata), rownames(nodeat_df))
            rdata <- cbind(rdata, nodeat_df[idx_cbind, , drop = FALSE])
        }
        
        return(rdata)
    }
)

#' @rdname GraphExperiment-methods
#' @importFrom SingleCellExperiment rowData
#' @export
setMethod(
    "graphs", "GraphExperiment", 
    function(x) {
        glist <- x@graphs
        if(length(glist) >0) {
            rdata <- rowData(x)
            glist <- SimpleList(lapply(glist, function(x) {
                .rowdata2nat(rdata, x)
            }))
        }
        return(glist)
    }
)


#' @rdname GraphExperiment-methods
#' @export
setMethod(
    "graph", c("GraphExperiment", "missing"),
    function(x, i) {
        
        glist <- x@graphs
        if(length(glist) == 0) {
            stop("The 'graphs' slot is empty.")
        }
        .rowdata2nat(rowData(x), glist[[1]])
    }
)

#' @rdname GraphExperiment-methods
#' @importFrom igraph set_vertex_attr
#' @export
setMethod(
    "graph", "GraphExperiment", 
    function(x, i) {
        
        tryCatch({
            .rowdata2nat(rowData(x), x@graphs[[i]])
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
