


#' @importFrom methods callNextMethod
#' @importFrom S4Vectors coolcat
setMethod(
    "show", "GraphExperiment",
    function(object) {
        callNextMethod()
        coolcat("graphs(%d): %s\n", graphNames(object))
    }
)