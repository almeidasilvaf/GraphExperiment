


#' @importFrom methods callNextMethod
#' @importFrom S4Vectors coolcat
setMethod(
    "show", "GraphExperiment",
    function(object) {
        callNextMethod()
        coolcat("rowGraphs(%d): %s\n", rowGraphNames(object))
        coolcat("colGraphs(%d): %s\n", colGraphNames(object))
    }
)