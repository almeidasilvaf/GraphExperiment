# `GraphExperiment` S4 class

The `GraphExperiment` class was designed to represent rectangular,
quantitative data (e.g., from transcriptomics, proteomics, metabolomics)
along with graphs showing how features (e.g., genes, proteins,
compounds) and observations (e.g., samples, cells, species, spots)
interact with each other. It extends `SingleCellExperiment` by providing
users with additional slots where row/column graphs can be stored.

## Usage

``` r
GraphExperiment(..., rowGraphs = list(), colGraphs = list())
```

## Arguments

- ...:

  Arguments passed to the `SingleCellExperiment` constructor function.

- rowGraphs:

  A list of `igraph` objects with one or multiple graphs representing
  how features relate to each other. Node names (i.e., `V(graph)$name`)
  must match rownames of assays.

- colGraphs:

  A list of `igraph` objects with one or multiple graphs representing
  how observations relate to each other. Node names (i.e.,
  `V(graph)$name`) must match colnames of assays.

## Value

A `GraphExperiment` object.

## Details

Like `SingleCellExperiment`, the `GraphExperiment` S4 class stores
quantitative data with associated metadata (i.e., `rowData` and
`colData`) along with embeddings from dimensionality reduction
techniques. However, it provides users with additional containers for
`igraph` objects containing graphs describing how features and/or
observations interact with each other. Graphs for features are stored in
a `rowGraphs` slot, and graphs for observations are stored in a
`colGraphs` slot. Both slots can hold one or multiple `igraph` objects
with some sort of network representation of the features in `rownames`
or observations in `colnames`. Example graphs for features can be
coexpression networks, regulatory networks, or co-abundance networks.
Example graphs for observations can be cell-cell distances,
sample-sample distances, or species networks representing genealogies
(as an alternative to phylogenies).

Importantly, node names in each row graph must match `rownames` of the
assays, and node names in each column graph must match `colnames` of the
assays. Besides, subsetting methods simultaneously subset `assays`,
`rowData`, `rowGraphs`, `colData`, and `colGraphs`.

Besides the constructor function (`GraphExperiment()`), a
`GraphExperiment` object can also be created by coercing from a
`SummarizedExperiment` or `SingleCellExperiment` object.

## Examples

``` r
# Example 1: from constructor function ----
## Simulate a matrix with 200 genes and 100 cells
gene_ids <- paste0("gene", seq_len(200))
cell_ids <- paste0("cell", seq_len(100))
mat <- matrix(rpois(20000, 5), ncol = 100, dimnames = list(gene_ids, cell_ids))

## Create a rowGraph from correlations (`igraph` object)
g <- graph_from_adjacency_matrix(cor(t(mat)), weighted = TRUE)

## Create a colGraph, also from correlations (but see scran::buildSNNGraph)
g2 <- graph_from_adjacency_matrix(cor(mat), weighted = TRUE)

## Construct `GraphExperiment` object
ge <- GraphExperiment(
    assays = list(counts = mat),
    rowGraphs = list(cor = g),
    colGraphs = list(cellcor = g2)
)
ge
#> class: GraphExperiment 
#> dim: 200 100 
#> metadata(0):
#> assays(1): counts
#> rownames(200): gene1 gene2 ... gene199 gene200
#> rowData names(0):
#> colnames(100): cell1 cell2 ... cell99 cell100
#> colData names(0):
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
#> rowGraphs(1): cor
#> colGraphs(1): cellcor

# Example 2: From `SingleCellExperiment` object ----
sce <- SingleCellExperiment(assays = list(counts = mat))
ge <- as(sce, "GraphExperiment")
ge
#> class: GraphExperiment 
#> dim: 200 100 
#> metadata(0):
#> assays(1): counts
#> rownames(200): gene1 gene2 ... gene199 gene200
#> rowData names(0):
#> colnames(100): cell1 cell2 ... cell99 cell100
#> colData names(0):
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
#> rowGraphs(0):
#> colGraphs(0):
```
