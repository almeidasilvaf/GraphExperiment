# Methods for `GraphExperiment` objects

The `GraphExperiment` class provides users with methods to get and set
graphs (`igraph` objects) representing how features and observations of
`SingleCellExperiment` objects relate to each other.

## Usage

``` r
# S4 method for class 'GraphExperiment'
rowData(x, use.names = TRUE, ...)

# S4 method for class 'GraphExperiment'
colData(x, use.names = TRUE, ...)

# S4 method for class 'GraphExperiment'
rowGraphs(x)

# S4 method for class 'GraphExperiment'
colGraphs(x)

# S4 method for class 'GraphExperiment,missing'
rowGraph(x, i)

# S4 method for class 'GraphExperiment,ANY'
rowGraph(x, i)

# S4 method for class 'GraphExperiment,missing'
colGraph(x, i)

# S4 method for class 'GraphExperiment,ANY'
colGraph(x, i)

# S4 method for class 'GraphExperiment'
rowGraphNames(x)

# S4 method for class 'GraphExperiment'
colGraphNames(x)

# S4 method for class 'GraphExperiment'
rowGraphs(x) <- value

# S4 method for class 'GraphExperiment'
colGraphs(x) <- value

# S4 method for class 'GraphExperiment'
rowGraph(x, i) <- value

# S4 method for class 'GraphExperiment'
colGraph(x, i) <- value

# S4 method for class 'GraphExperiment,character'
rowGraphNames(x) <- value

# S4 method for class 'GraphExperiment,character'
colGraphNames(x) <- value
```

## Arguments

- x:

  A `GraphExperiment` object.

- use.names:

  Passed to the `colData` method of `SingleCellExperiment`. Default:
  TRUE.

- ...:

  Ignored.

- i:

  List element (numeric for index, character for name) of the element to
  access or replace.

- value:

  Replacement value for replacement methods.

## Value

Return values depend on the method. See details and examples.

## rowGraphs/colGraphs and rowGraph/colGraph methods

- `rowGraphs(x)` and `colGraphs(x)`: :

  Getter for a `SimpleList` of `igraph` objects representing rows and
  columns, respectively.

- `rowGraphs(x) <- value` and `colGraphs(x) <- value`: :

  Setter for a SimpleList or list (coerced to SimpleList) of `igraph`
  objects representing rows and columns, respectively.

- `rowGraph(x, i)` and `colGraph(x, i)`: :

  Getter for an `igraph` object containing graph `i` from the list
  stored in `rowGraphs` and `colGraphs`, respectively.

- `rowGraph(x, i) <- value` and `colGraph(x, i) <- value`: :

  Setter for an `igraph` object to be stored in element `i` of
  `rowGraphs` and `colGraphs`.

## rowGraphNames and colGraphNames methods

- `rowGraphNames(x)` and `colGraphNames(x)`: :

  Getter to extract names of graphs in `rowGraphs` and `colGraphs`.

- `rowGraphNames(x) <- value` and `colGraphNames(x) <- value`: :

  Setter to assign new names to the graphs stored in `rowGraphs` and
  `colGraphs`.

## rowData method

- `rowData(x)`: :

  Getter to extract rowData (as in `SingleCellExperiment` objects), but
  with node attributes of graphs included.

## colData method

- `colData(x)`: :

  Getter to extract colData (as in `SingleCellExperiment` objects), but
  with node attributes of graphs included.

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

# Extract graph names
rowGraphNames(ge)
#> [1] "cor"
colGraphNames(ge)
#> [1] "cellcor"


# Extract graphs
rowGraphs(ge)
#> List of length 1
#> names(1): cor
rowGraph(ge, "cor")
#> IGRAPH ce6caa6 DNW- 200 39990 -- 
#> + attr: name (v/c), degree (v/n), pathway (v/c), coding (v/l), weight
#> | (e/n)
#> + edges from ce6caa6 (vertex names):
#>  [1] gene1->gene1  gene1->gene2  gene1->gene3  gene1->gene4  gene1->gene5 
#>  [6] gene1->gene6  gene1->gene7  gene1->gene8  gene1->gene9  gene1->gene10
#> [11] gene1->gene11 gene1->gene12 gene1->gene13 gene1->gene14 gene1->gene15
#> [16] gene1->gene16 gene1->gene17 gene1->gene18 gene1->gene19 gene1->gene20
#> [21] gene1->gene21 gene1->gene22 gene1->gene23 gene1->gene24 gene1->gene25
#> [26] gene1->gene26 gene1->gene27 gene1->gene28 gene1->gene29 gene1->gene30
#> [31] gene1->gene31 gene1->gene32 gene1->gene33 gene1->gene34 gene1->gene35
#> + ... omitted several edges

colGraphs(ge)
#> List of length 1
#> names(1): cellcor
colGraph(ge, "cellcor")
#> IGRAPH 0460ba5 DNW- 100 9998 -- 
#> + attr: name (v/c), celltype (v/c), weight (e/n)
#> + edges from 0460ba5 (vertex names):
#>  [1] cell1->cell1  cell1->cell2  cell1->cell3  cell1->cell4  cell1->cell5 
#>  [6] cell1->cell6  cell1->cell7  cell1->cell8  cell1->cell9  cell1->cell10
#> [11] cell1->cell11 cell1->cell12 cell1->cell13 cell1->cell14 cell1->cell15
#> [16] cell1->cell16 cell1->cell17 cell1->cell18 cell1->cell19 cell1->cell20
#> [21] cell1->cell21 cell1->cell22 cell1->cell23 cell1->cell24 cell1->cell25
#> [26] cell1->cell26 cell1->cell27 cell1->cell28 cell1->cell29 cell1->cell30
#> [31] cell1->cell31 cell1->cell32 cell1->cell33 cell1->cell34 cell1->cell35
#> [36] cell1->cell36 cell1->cell37 cell1->cell38 cell1->cell39 cell1->cell40
#> + ... omitted several edges


# Add a new graph
rowGraph(ge, "newcor") <- g
ge
#> class: GraphExperiment 
#> dim: 200 100 
#> metadata(0):
#> assays(1): counts
#> rownames(200): gene1 gene2 ... gene199 gene200
#> rowData names(4): pathway coding cor__degree newcor__degree
#> colnames(100): cell1 cell2 ... cell99 cell100
#> colData names(1): celltype
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
#> rowGraphs(2): cor newcor
#> colGraphs(1): cellcor

# Add a list of graphs
colGraphs(ge) <- list(cellcor = g2, new_cellcor = g2)
ge
#> class: GraphExperiment 
#> dim: 200 100 
#> metadata(0):
#> assays(1): counts
#> rownames(200): gene1 gene2 ... gene199 gene200
#> rowData names(4): pathway coding cor__degree newcor__degree
#> colnames(100): cell1 cell2 ... cell99 cell100
#> colData names(1): celltype
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
#> rowGraphs(2): cor newcor
#> colGraphs(2): cellcor new_cellcor

# Replace graph names
rowGraphNames(ge) <- c("network1", "network2")
ge
#> class: GraphExperiment 
#> dim: 200 100 
#> metadata(0):
#> assays(1): counts
#> rownames(200): gene1 gene2 ... gene199 gene200
#> rowData names(4): pathway coding network1__degree network2__degree
#> colnames(100): cell1 cell2 ... cell99 cell100
#> colData names(1): celltype
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
#> rowGraphs(2): network1 network2
#> colGraphs(2): cellcor new_cellcor

# Access rowData (note: rowData + node attributes combined)
rowData(ge)
#> DataFrame with 200 rows and 4 columns
#>             pathway    coding network1__degree network2__degree
#>         <character> <logical>        <numeric>        <numeric>
#> gene1            P1      TRUE         3.015413         3.015413
#> gene2            P1      TRUE         0.751915         0.751915
#> gene3            P2     FALSE         3.571644         3.571644
#> gene4            P1      TRUE         6.074124         6.074124
#> gene5            P1      TRUE         1.342499         1.342499
#> ...             ...       ...              ...              ...
#> gene196          P2      TRUE        -0.387201        -0.387201
#> gene197          P1      TRUE         0.789255         0.789255
#> gene198          P1      TRUE         6.323828         6.323828
#> gene199          P2      TRUE         5.638062         5.638062
#> gene200          P2     FALSE        -2.567450        -2.567450

# Access colData (note: colData + node attributes combined)
colData(ge)
#> DataFrame with 100 rows and 1 column
#>            celltype
#>         <character>
#> cell1           ct2
#> cell2           ct1
#> cell3           ct2
#> cell4           ct2
#> cell5           ct1
#> ...             ...
#> cell96          ct1
#> cell97          ct2
#> cell98          ct1
#> cell99          ct1
#> cell100         ct2
```
