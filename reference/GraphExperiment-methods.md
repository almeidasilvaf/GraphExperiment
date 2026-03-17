# Methods for `GraphExperiment` objects

The `GraphExperiment` class provides users with methods to get and set
graphs (`igraph` objects) representing how features of
`SingleCellExperiment` objects relate to each other.

## Usage

``` r
# S4 method for class 'GraphExperiment'
rowData(x, use.names = TRUE, ...)

# S4 method for class 'GraphExperiment'
graphs(x)

# S4 method for class 'GraphExperiment,missing'
graph(x, i)

# S4 method for class 'GraphExperiment,ANY'
graph(x, i)

# S4 method for class 'GraphExperiment'
graphNames(x)

# S4 method for class 'GraphExperiment'
graphs(x) <- value

# S4 method for class 'GraphExperiment'
graph(x, i) <- value

# S4 method for class 'GraphExperiment,character'
graphNames(x) <- value
```

## Arguments

- x:

  A `GraphExperiment` object.

- use.names:

  Passed to the `rowData` method of `SingleCellExperiment`. Default:
  TRUE.

- ...:

  Ignored.

- i:

  List element (numeric for index, character for name) of the element to
  access or replace.

- value:

  Replacement value for replacement methods

## Value

Return values depend on the method. See details and examples.

## graphs and graph methods

- `graphs(x)`: :

  Getter for a `SimpleList` of `igraph` objects.

- `graphs(x) <- value`: :

  Setter for a SimpleList or list (coerced to SimpleList) of `igraph`
  objects.

- `graph(x, i)`: :

  Getter for an `igraph` object containing graph `i` from the list
  stored in `graphs`.

- `graph(x, i) <- value`: :

  Setter for an `igraph` object to be stored in element `i` of the list
  `graphs`

## graphNames methods

- `graphNames(x)`: :

  Getter to extract names of graphs in the `graphs` slot.

- `graphNames(x) <- value`: :

  Setter to assign new names to the graphs stored in the `graphs` slot.

## Examples

``` r
# Simulate elements of a GraphExperiment object
## Assays
gene_ids <- paste0("gene", seq_len(200))
cell_ids <- paste0("cell", seq_len(100))
mat <- matrix(rpois(20000, 5), ncol = 100, dimnames = list(gene_ids, cell_ids))

## Graph (with node attributes)
g <- graph_from_adjacency_matrix(cor(t(mat)), weighted = TRUE)
V(g)$degree <- igraph::strength(g)

## rowData
rdata <- data.frame(
    row.names = gene_ids, 
    pathway = sample(c("P1", "P2"), size = length(gene_ids), replace = TRUE),
    coding = sample(c(TRUE, FALSE), size = length(gene_ids), replace = TRUE)
)


# Create a GraphExperiment object
ge <- GraphExperiment(
    assays = list(counts = mat), 
    rowData = rdata,
    graphs = list(cor = g)
)
ge
#> class: GraphExperiment 
#> dim: 200 100 
#> metadata(0):
#> assays(1): counts
#> rownames(200): gene1 gene2 ... gene199 gene200
#> rowData names(3): pathway coding cor__degree
#> colnames(100): cell1 cell2 ... cell99 cell100
#> colData names(0):
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
#> graphs(1): cor

# Extract graph names
graphNames(ge)
#> [1] "cor"
ge
#> class: GraphExperiment 
#> dim: 200 100 
#> metadata(0):
#> assays(1): counts
#> rownames(200): gene1 gene2 ... gene199 gene200
#> rowData names(3): pathway coding cor__degree
#> colnames(100): cell1 cell2 ... cell99 cell100
#> colData names(0):
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
#> graphs(1): cor

# Extract graphs
graphs(ge)
#> List of length 1
#> names(1): cor
graph(ge, "cor")
#> IGRAPH e2e818c DNW- 200 39990 -- 
#> + attr: name (v/c), degree (v/n), pathway (v/c), coding (v/l), weight
#> | (e/n)
#> + edges from e2e818c (vertex names):
#>  [1] gene1->gene1  gene1->gene2  gene1->gene3  gene1->gene4  gene1->gene5 
#>  [6] gene1->gene6  gene1->gene7  gene1->gene8  gene1->gene9  gene1->gene10
#> [11] gene1->gene11 gene1->gene12 gene1->gene13 gene1->gene14 gene1->gene15
#> [16] gene1->gene16 gene1->gene17 gene1->gene18 gene1->gene19 gene1->gene20
#> [21] gene1->gene21 gene1->gene22 gene1->gene23 gene1->gene24 gene1->gene25
#> [26] gene1->gene26 gene1->gene27 gene1->gene28 gene1->gene29 gene1->gene30
#> [31] gene1->gene31 gene1->gene32 gene1->gene33 gene1->gene34 gene1->gene35
#> + ... omitted several edges

# Add a new graph
graph(ge, "newcor") <- g
ge
#> class: GraphExperiment 
#> dim: 200 100 
#> metadata(0):
#> assays(1): counts
#> rownames(200): gene1 gene2 ... gene199 gene200
#> rowData names(6): pathway coding ... cor__coding newcor__degree
#> colnames(100): cell1 cell2 ... cell99 cell100
#> colData names(0):
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
#> graphs(2): cor newcor

# Add a list of graphs
graphs(ge) <- list(cor = g, newcor = g)
ge
#> class: GraphExperiment 
#> dim: 200 100 
#> metadata(0):
#> assays(1): counts
#> rownames(200): gene1 gene2 ... gene199 gene200
#> rowData names(4): pathway coding cor__degree newcor__degree
#> colnames(100): cell1 cell2 ... cell99 cell100
#> colData names(0):
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
#> graphs(2): cor newcor

# Replace graph names
graphNames(ge) <- c("network1", "network2")
ge
#> class: GraphExperiment 
#> dim: 200 100 
#> metadata(0):
#> assays(1): counts
#> rownames(200): gene1 gene2 ... gene199 gene200
#> rowData names(8): pathway coding ... network2__pathway network2__coding
#> colnames(100): cell1 cell2 ... cell99 cell100
#> colData names(0):
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
#> graphs(2): network1 network2

# Access rowData (note: rowData + node attributes combined)
rowData(ge)
#> DataFrame with 200 rows and 8 columns
#>             pathway    coding network1__degree network1__pathway
#>         <character> <logical>        <numeric>       <character>
#> gene1            P1      TRUE         3.015413                P1
#> gene2            P1      TRUE         0.751915                P1
#> gene3            P2     FALSE         3.571644                P2
#> gene4            P1      TRUE         6.074124                P1
#> gene5            P1      TRUE         1.342499                P1
#> ...             ...       ...              ...               ...
#> gene196          P2      TRUE        -0.387201                P2
#> gene197          P1      TRUE         0.789255                P1
#> gene198          P1      TRUE         6.323828                P1
#> gene199          P2      TRUE         5.638062                P2
#> gene200          P2     FALSE        -2.567450                P2
#>         network1__coding network2__degree network2__pathway network2__coding
#>                <logical>        <numeric>       <character>        <logical>
#> gene1               TRUE         3.015413                P1             TRUE
#> gene2               TRUE         0.751915                P1             TRUE
#> gene3              FALSE         3.571644                P2            FALSE
#> gene4               TRUE         6.074124                P1             TRUE
#> gene5               TRUE         1.342499                P1             TRUE
#> ...                  ...              ...               ...              ...
#> gene196             TRUE        -0.387201                P2             TRUE
#> gene197             TRUE         0.789255                P1             TRUE
#> gene198             TRUE         6.323828                P1             TRUE
#> gene199             TRUE         5.638062                P2             TRUE
#> gene200            FALSE        -2.567450                P2            FALSE
```
