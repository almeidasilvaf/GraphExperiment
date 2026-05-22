
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
        rowGraphs = "list_or_SimpleList",
        colGraphs = "list_or_SimpleList"
    ),
    prototype = prototype(rowGraphs = list(), colGraphs = list())
)