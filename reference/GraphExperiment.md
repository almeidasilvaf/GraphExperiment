# `GraphExperiment` S4 class

The `GraphExperiment` class was designed to represent rectangular,
quantitative data (e.g., from transcriptomics, proteomics, metabolomics)
along with graphs showing how features (e.g., genes, proteins,
compounds) interact with each other. It extends `SingleCellExperiment`
by providing users with an additional slot where graphs can be stored.

## Usage

``` r
GraphExperiment(..., graphs = list())
```

## Arguments

- ...:

  Arguments passed to the `SingleCellExperiment` constructor function.

- graphs:

  A list of `igraph` objects with one or multiple graphs. Node names
  (i.e., `V(graph)$name`) must match rownames of assays.

## Value

A `GraphExperiment` object.

## Details

Like `SingleCellExperiment`, the `GraphExperiment` S4 class stores
quantitative data with associated metadata (i.e., `rowData` and
`colData`) along with embeddings from dimensionality reduction
techniques. However, it provides users with an additional container for
`igraph` objects containing graphs describing how features interact with
each other. Graphs are stored in a `graphs` slot, which can hold one or
multiple `igraph` objects with some sort of network representation of
the features in `rownames`. Example graphs can be coexpression networks,
regulatory networks, or co-abundance networks.

Importantly, node names in each graph must match `rownames` of the
assays, and subsetting methods simultaneously subset `assays`,
`rowData`, and `graphs`.

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

## Create a graph from correlations (`igraph` object)
g <- graph_from_adjacency_matrix(cor(t(mat)), weighted = TRUE)

## Construct `GraphExperiment` object
ge <- GraphExperiment(
    assays = list(counts = mat),
    graphs = list(cor = g)
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
#> graphs(1): cor

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
#> graphs(0):
```
