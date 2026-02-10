
## Getters ---------------------------------------------------------------------

#' @export
setGeneric("graphs", function(x, ...) standardGeneric("graphs"))

#' @export
setGeneric(
    "graph", signature = c("x", "i"),
    function(x, i, ...) standardGeneric("graph")
)

#' @export
setGeneric("graphNames", function(x, ...) standardGeneric("graphNames"))


## Setters ---------------------------------------------------------------------

#' @export
setGeneric("graphs<-", function(x, ..., value) standardGeneric("graphs<-"))

#' @export
setGeneric(
    "graph<-", signature = c("x", "i"),
    function(x, i, ..., value) standardGeneric("graph<-")
)

#' @export
setGeneric("graphNames<-", function(x, ..., value) standardGeneric("graphNames<-"))
