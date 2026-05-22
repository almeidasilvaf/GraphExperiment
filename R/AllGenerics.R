
## Getters ---------------------------------------------------------------------

#' @export
setGeneric("rowGraphs", function(x, ...) standardGeneric("rowGraphs"))

#' @export
setGeneric("colGraphs", function(x, ...) standardGeneric("colGraphs"))




#' @export
setGeneric(
    "rowGraph", signature = c("x", "i"),
    function(x, i, ...) standardGeneric("rowGraph")
)

#' @export
setGeneric(
    "colGraph", signature = c("x", "i"),
    function(x, i, ...) standardGeneric("colGraph")
)




#' @export
setGeneric("rowGraphNames", function(x, ...) standardGeneric("rowGraphNames"))

#' @export
setGeneric("colGraphNames", function(x, ...) standardGeneric("colGraphNames"))


## Setters ---------------------------------------------------------------------

#' @export
setGeneric("rowGraphs<-", function(x, ..., value) standardGeneric("rowGraphs<-"))

#' @export
setGeneric("colGraphs<-", function(x, ..., value) standardGeneric("colGraphs<-"))




#' @export
setGeneric(
    "rowGraph<-", signature = c("x", "i"),
    function(x, i, ..., value) standardGeneric("rowGraph<-")
)

#' @export
setGeneric(
    "colGraph<-", signature = c("x", "i"),
    function(x, i, ..., value) standardGeneric("colGraph<-")
)




#' @export
setGeneric("rowGraphNames<-", function(x, ..., value) standardGeneric("rowGraphNames<-"))

#' @export
setGeneric("colGraphNames<-", function(x, ..., value) standardGeneric("colGraphNames<-"))

