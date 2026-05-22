# Subsetting `GraphExperiment` objects

The subsetting method for [`GraphExperiment`](GraphExperiment.md)
objects ensures that nodes from `igraph` objects are filtered to match
rows and columns with the remainder of the object.

## Arguments

- x:

  A [`GraphExperiment`](GraphExperiment.md) object.

- i:

  Numeric, row indices for subsetting.

- j:

  Numeric, column indices for subsetting.

## Value

a [`GraphExperiment`](GraphExperiment.md) object.

## subset

- `[`::

  subsetting method

## Examples

``` r
# Simulate elements of a GraphExperiment object
## Assays
gene_ids <- paste0("gene", seq_len(200))
cell_ids <- paste0("cell", seq_len(100))
mat <- matrix(rpois(20000, 5), ncol = 100, dimnames = list(gene_ids, cell_ids))

## rowGraph (with node attributes)
g <- graph_from_adjacency_matrix(cor(t(mat)), weighted = TRUE)
V(g)$degree <- igraph::strength(g)

## colGraph
g2 <- graph_from_adjacency_matrix(cor(mat), weighted = TRUE)

## rowData
rdata <- data.frame(
    row.names = gene_ids, 
    pathway = sample(c("P1", "P2"), size = length(gene_ids), replace = TRUE),
    coding = sample(c(TRUE, FALSE), size = length(gene_ids), replace = TRUE)
)

## colData
cdata <- data.frame(
    row.names = cell_ids, 
    celltype = sample(c("ct1", "ct2"), size = length(cell_ids), replace = TRUE)
)


# Create a GraphExperiment object
ge <- GraphExperiment(
    assays = list(counts = mat), 
    rowData = rdata,
    colData = cdata,
    rowGraphs = list(cor = g),
    colGraphs = list(cellcor = g2)
)
ge
#> class: GraphExperiment 
#> dim: 200 100 
#> metadata(0):
#> assays(1): counts
#> rownames(200): gene1 gene2 ... gene199 gene200
#> rowData names(3): pathway coding cor__degree
#> colnames(100): cell1 cell2 ... cell99 cell100
#> colData names(1): celltype
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
#> rowGraphs(1): cor
#> colGraphs(1): cellcor

# Subset object
ge[1:5, 1:5]
#> class: GraphExperiment 
#> dim: 5 5 
#> metadata(0):
#> assays(1): counts
#> rownames(5): gene1 gene2 gene3 gene4 gene5
#> rowData names(3): pathway coding cor__degree
#> colnames(5): cell1 cell2 cell3 cell4 cell5
#> colData names(1): celltype
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
#> rowGraphs(1): cor
#> colGraphs(1): cellcor
```
