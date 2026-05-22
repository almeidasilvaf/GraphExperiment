
#' @importFrom igraph V
.nodes_validity_rowgraphs <- function(object) {
    
    graph_list <- object@rowGraphs
    if(length(graph_list) == 0) return(NULL)
    rn <- rownames(object)
    
    msg <- NULL
    msgs <- unlist(lapply(seq_along(graph_list), function(i) {
        g <- graph_list[[i]]
        graph_name <- names(graph_list)[i]
        
        node_names <- V(g)$name
        
        # Check for features in rownames but not in graph
        missing_from_graph <- setdiff(rn, node_names)
        if(length(missing_from_graph) > 0) {
            msg <- c(msg, sprintf(
                "%d feature(s) in 'rownames' are missing from graph '%s'.",
                length(missing_from_graph), graph_name
            ))
        }
        
        # Check for nodes in graph but not in rownames
        missing_from_rownames <- setdiff(node_names, rn)
        if(length(missing_from_rownames) > 0) {
            msg <- c(msg, sprintf(
                "%d node(s) in graph '%s' are missing from 'rownames'.",
                length(missing_from_rownames), graph_name
            ))
        }
        
        return(msg)
    }))
    
    return(msgs)
}

.nodes_validity_colgraphs <- function(object) {
    
    graph_list <- object@colGraphs
    if(length(graph_list) == 0) return(NULL)
    cn <- colnames(object)
    
    msg <- NULL
    msgs <- unlist(lapply(seq_along(graph_list), function(i) {
        g <- graph_list[[i]]
        graph_name <- names(graph_list)[i]
        
        node_names <- V(g)$name
        
        # Check for observations in colnames but not in graph
        missing_from_graph <- setdiff(cn, node_names)
        if(length(missing_from_graph) > 0) {
            msg <- c(msg, sprintf(
                "%d observations(s) in 'colnames' are missing from graph '%s'.",
                length(missing_from_graph), graph_name
            ))
        }
        
        # Check for observations in graph but not in colnames
        missing_from_rownames <- setdiff(node_names, cn)
        if(length(missing_from_rownames) > 0) {
            msg <- c(msg, sprintf(
                "%d node(s) in graph '%s' are missing from 'colnames'.",
                length(missing_from_rownames), graph_name
            ))
        }
        
        return(msg)
    }))
    
    return(msgs)
}


#' @importFrom methods is
.graphs_validity <- function(object) {
    
    rg <- object@rowGraphs
    cg <- object@colGraphs
    
    # Check that each element is an igraph object
    if(length(rg) > 0) {
        msgs <- unlist(lapply(seq_along(rg), function(i) {
            graph_name <- names(rg)[i]
            if(!is(rg[[i]], "igraph")) {
                return(sprintf("Graph '%s' in rowGraphs is not an igraph object", graph_name))
            }
        }))
    }
    
    if(length(cg) > 0) {
        msgs <- unlist(lapply(seq_along(cg), function(i) {
            graph_name <- names(cg)[i]
            if(!is(cg[[i]], "igraph")) {
                return(sprintf("Graph '%s' in colGraphs is not an igraph object", graph_name))
            }
        }))
    }
    
    return(NULL)
}


.ge_validity <- function(object) {
    
    msg <- .graphs_validity(object)
    if(!is.null(msg)) return(msg)
    
    msg <- .nodes_validity_rowgraphs(object)
    if(!is.null(msg)) return(msg)
    
    msg <- .nodes_validity_colgraphs(object)
    if(!is.null(msg)) return(msg)
    
    return(TRUE)
}


#' @importFrom S4Vectors setValidity2
setValidity2("GraphExperiment", .ge_validity)
