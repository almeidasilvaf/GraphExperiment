
#' Coercion methods for `GraphExperiment` objects
#' 
#' The \code{GraphExperiment} class inherits from \code{SingleCellExperiment},
#' which in turn inherits from \code{(Ranged)SummarizedExperiment}. To ensure
#' \code{GraphExperiment} easily interoperates with these classes, we provide
#' users with traditional \code{as()} coercion methods.
#'
#' @param from An object of class \code{SingleCellExperiment}, 
#' \code{SummarizedExperiment}, or \code{RangedSummarizedExperiment} to be 
#' coerced to \code{GraphExperiment}.
#' 
#' @return A \code{GraphExperiment} object.
#' 
#' @name GraphExperiment-coerce
#' @aliases 
#' coerce,SingleCellExperiment,GraphExperiment-method
#' coerce,SummarizedExperiment,GraphExperiment-method
#' coerce,RangedSummarizedExperiment,GraphExperiment-method
#' 
#' @examples 
#' # Simulate a count matrix
#' gene_ids <- paste0("gene", seq_len(200))
#' cell_ids <- paste0("cell", seq_len(100))
#' mat <- matrix(rpois(20000, 5), ncol = 100, dimnames = list(gene_ids, cell_ids))
#' 
#' # Coerce from `SummarizedExperiment`
#' se <- SummarizedExperiment(assays = list(mat))
#' as(se, "GraphExperiment")
#' 
#' # Coerce from `RangedSummarizedExperiment`
#' rse <- as(se, "RangedSummarizedExperiment")
#' as(rse, "GraphExperiment")
#' 
#' # Coerce from `SingleCellExperiment`
#' sce <- as(se, "SingleCellExperiment")
#' as(sce, "GraphExperiment")
NULL


#' @export
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
setAs(
    "SingleCellExperiment", "GraphExperiment", 
    function(from) {
        new("GraphExperiment", from, graphs = SimpleList())
    }
)

#' @export
#' @importClassesFrom SummarizedExperiment SummarizedExperiment 
setAs(
    "SummarizedExperiment", "GraphExperiment", 
    function(from) {
        new("GraphExperiment", as(from, "SingleCellExperiment"), graphs = SimpleList())
    }
)

#' @export
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment 
setAs(
    "RangedSummarizedExperiment", "GraphExperiment", 
    function(from) {
        new("GraphExperiment", as(from, "SingleCellExperiment"), graphs = SimpleList())
    }
)
