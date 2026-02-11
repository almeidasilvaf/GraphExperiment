
#' @importFrom igraph V
.nodes_validity <- function(object) {
    
    graph_list <- graphs(object)
    if(length(graph_list) == 0) return(NULL)
    rn <- rownames(object)
    
    msg <- NULL
    msgs <- unlist(lapply(seq_along(graph_list), function(i) {
        g <- graph_list[[i]]
        graph_name <- names(graph_list)[i]
        if(is.null(graph_name)) graph_name <- as.character(i)
        
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

#' @importFrom methods is
.graphs_validity <- function(object) {
    
    g <- graphs(object)
    
    # Check that graphs slot has correct type
    if(!is.null(g) && !is(g, "list") && !is(g, "SimpleList")) {
        return("'graphs' must be a list or SimpleList.")
    }
    
    # Check that each element is an igraph object
    if(length(g) > 0) {
        msgs <- unlist(lapply(seq_along(g), function(i) {
            graph_name <- names(g)[i]
            if(is.null(graph_name)) graph_name <- as.character(i)
            if(!is(g[[i]], "igraph")) {
                return(sprintf("Graph '%s' is not an igraph object", graph_name))
            }
        }))
    }
    
    return(NULL)
}

.ge_validity <- function(object) {
    
    msg <- .graphs_validity(object)
    if(!is.null(msg)) return(msg)
    
    msg <- .nodes_validity(object)
    if(!is.null(msg)) return(msg)
    
    return(TRUE)
}


#' @importFrom S4Vectors setValidity2
setValidity2("GraphExperiment", .ge_validity)
