
setClassUnion("list_or_SimpleList", c("list", "SimpleList", "NULL"))

#' `GraphExperiment` S4 class
#' 
#' @export
#' @rdname GraphExperiment
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importClassesFrom S4Vectors SimpleList
setClass(
    "GraphExperiment",
    contains = "SingleCellExperiment",
    slots = c(
        graphs = "list_or_SimpleList"
    ),
    prototype = prototype(graph = list())
)