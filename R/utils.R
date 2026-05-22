
#' Add rowData or colData variables to node attributes of `igraph` objects
#' 
#' @importFrom igraph set_vertex_attr vertex_attr_names
#' @noRd
.metadata2nat <- function(metadata, graph) {
    
    # Remove rowData vars extracted from graph node attributes
    idx_remove <- grepl("__", names(metadata))
    if(sum(idx_remove) >0) { metadata <- metadata[, !idx_remove, drop = FALSE] }
    
    # Add missing row/coldata vars to node attributes
    fgraph <- graph
    if(length(metadata) >0) {
        nat_list <- as.list(metadata)
        
        if(length(nat_list) >0) {
            idx <- match(rownames(metadata), V(graph)$name)
            for(rvar in names(nat_list)) {
                fgraph <- set_vertex_attr(
                    fgraph, name = rvar, index = idx, value = nat_list[[rvar]]
                )
            }
        }
    }
    
    return(fgraph)
}


#' Retrieve rowData/colData that combines node attributes of `igraph` object
#' 
#' @importFrom igraph as_data_frame
#' @noRd
.get_metadata <- function(graph_list, metadata) {
    
    glist <- graph_list
    gnames <- names(glist)
    
    # Get a data frame of node attributes for all graphs (if any)
    nodeat_df <- NULL
    if(length(glist) > 0) {
        nodeat_df <- Reduce(cbind, lapply(seq_along(glist), function(n) {
            df <- igraph::as_data_frame(glist[[n]], what = "vertices")

            remove_idx <- which(names(df) %in% colnames(metadata)) 
            
            names(df) <- paste0(gnames[n], "__", names(df))
            df <- df[, -c(1, remove_idx), drop = FALSE]
            return(df)
        }))
    }
    
    # Extend row/colData to include node attributes
    if(length(nodeat_df) > 0) {
        idx_cbind <- match(rownames(metadata), rownames(nodeat_df))
        metadata <- cbind(metadata, nodeat_df[idx_cbind, , drop = FALSE])
    }
    
    return(metadata)
}


