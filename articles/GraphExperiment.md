# Introduction to the GraphExperiment class

## Introduction

Networks (or graphs) have become widely used data representations in
biology, as they can efficiently encode node-node interactions and
neighborhoods. In high-throughput, quantitative omics data (e.g.,
transcriptomics, proteomics, metabolomics, epigenomics, etc), widely
used feature-level network representations include gene coexpression,
protein-protein interaction, gene regulatory, and co-abundance networks,
and cell-cell communication, sample-sample similarity, and
species-species networks at the observation level. While data structures
to store quantitative data and associated metadata exist (e.g.,
`SummarizedExperiment`, `SingleCellExperiment`, `SpatialExperiment`,
etc), support for networks describing how features (i.e., assay rows)
and/or observations (i.e., assay columns) relate to each other is
currently missing. `GraphExperiment` is an S4 class that extends
`SingleCellExperiment` (Amezquita et al. 2020) to include additional
containers for networks associated with assay features (`rowGraphs`) and
assay observations (`colGraphs`).

Of note, trees are an alternative way of representing how assay
features/observations are related to each other. Users interested in
tree representations of assays rows/columns can use the
*[TreeSummarizedExperiment](https://bioconductor.org/packages/3.24/TreeSummarizedExperiment)*
package. Trees are essentially *a kind of graph* (i.e., all trees are
graphs, but not all graphs are trees). Here, we chose to use a more
general graph representation (namely `igraph` objects) to provide users
and developers with more flexibility.

## Installation

`GraphExperiment` can be installed from Bioconductor with the following
code:

``` r

if(!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install("GraphExperiment")
```

``` r

# Load package after installation
library(GraphExperiment)

set.seed(777) # for reproducibility
```

## Anatomy of a `GraphExperiment` object

Since the `GraphExperiment` class extends the `SingleCellExperiment`
class, all `SingleCellExperiment` slots are present in
`GraphExperiment`, including:

- `assays`: list of matrices with primary (e.g., counts) and transformed
  (e.g., log-normalized counts, TPM, etc) data, with features in rows
  and observations in columns.
- `colData`: a data frame with column (observation) metadata, such as
  cell type, sample ID, condition, batch ID, genotype, etc.
- `rowData`: a data frame with row (feature) metadata, such as gene ID,
  genomic coordinates, functional annotation, etc.
- `reducedDims`: list of data frames with reduced dimensions, such as
  PCA, t-SNE, and UMAP embeddings.

Compared to `SingleCellExperiment` objects, `GraphExperiment` provides
two additional containers: [^1]

- `rowGraphs`: list of `igraph` objects containing graphs representing
  feature-feature relationships, including (but optional) node and edge
  attributes.
- `colGraphs`: list of `igraph` objects containing graphs representing
  observation-observation relationships, including (but optional) node
  and edge attributes.

![The GraphExperiment class.](GraphExperiment.png)

The GraphExperiment class.

The `igraph` data class from the
*[igraph](https://CRAN.R-project.org/package=igraph)* package is the
standard data structure for graph representation in R. If you are
unfamiliar with `igraph` objects, you can learn more about it by reading
the *[igraph](https://CRAN.R-project.org/package=igraph)* vignettes.

## Building a `GraphExperiment` object

`GraphExperiment` objects can be created from scratch using the
constructor function
[`GraphExperiment()`](../reference/GraphExperiment.md). Below we will
simulate a scRNA-seq count matrix with some gene (row) and cell (column)
metadata, and create graphs based on gene-gene and cell-cell
correlations [^2].

``` r

# Simulate parts of a `GraphExperiment` object
## Assays
gene_ids <- paste0("gene", seq_len(200))
cell_ids <- paste0("cell", seq_len(100))
mat <- matrix(rpois(20000, 5), ncol = 100, dimnames = list(gene_ids, cell_ids))
mat[1:5, 1:5]
#>       cell1 cell2 cell3 cell4 cell5
#> gene1     6     7     5     5     6
#> gene2     5     3     6     3     5
#> gene3     4     5     9     7     2
#> gene4    12     4     5     9     3
#> gene5     6     6     5     2     8

## rowData
rdata <- data.frame(
    row.names = gene_ids,
    pathway = sample(c("P1", "P2"), size = length(gene_ids), replace = TRUE),
    coding = sample(c(TRUE, FALSE), size = length(gene_ids), replace = TRUE)
)
head(rdata)
#>       pathway coding
#> gene1      P1  FALSE
#> gene2      P1   TRUE
#> gene3      P2   TRUE
#> gene4      P1  FALSE
#> gene5      P1   TRUE
#> gene6      P1  FALSE

## colData
cdata <- data.frame(
    row.names = cell_ids, 
    cell_type = sample(c("ct1", "ct2"), size = length(cell_ids), replace = TRUE)
)
head(cdata)
#>       cell_type
#> cell1       ct2
#> cell2       ct2
#> cell3       ct2
#> cell4       ct1
#> cell5       ct2
#> cell6       ct1

## rowGraph (with node attribute `degree`)
rg <- graph_from_adjacency_matrix(
    cor(t(mat)), mode = "undirected", weighted = TRUE
)
rg <- set_vertex_attr(rg, "degree", value = strength(rg))
rg
#> IGRAPH 8c54443 UNW- 200 20096 -- 
#> + attr: name (v/c), degree (v/n), weight (e/n)
#> + edges from 8c54443 (vertex names):
#>  [1] gene1--gene1  gene1--gene2  gene1--gene3  gene1--gene4  gene1--gene5 
#>  [6] gene1--gene6  gene1--gene7  gene1--gene8  gene1--gene9  gene1--gene10
#> [11] gene1--gene11 gene1--gene12 gene1--gene13 gene1--gene14 gene1--gene15
#> [16] gene1--gene16 gene1--gene17 gene1--gene18 gene1--gene19 gene1--gene20
#> [21] gene1--gene21 gene1--gene22 gene1--gene23 gene1--gene24 gene1--gene25
#> [26] gene1--gene26 gene1--gene27 gene1--gene28 gene1--gene29 gene1--gene30
#> [31] gene1--gene31 gene1--gene32 gene1--gene33 gene1--gene34 gene1--gene35
#> [36] gene1--gene36 gene1--gene37 gene1--gene38 gene1--gene39 gene1--gene40
#> + ... omitted several edges

## colGraph
cg <- graph_from_adjacency_matrix(
    cor(mat), mode = "undirected", weighted = TRUE
)
```

To create a `GraphExperiment` object from the constructor function, you
would run:

``` r

# Create a `GraphExperiment` object
ge <- GraphExperiment(
    assays = list(counts = mat),
    rowData = rdata,
    colData = cdata,
    rowGraphs = list(gene_cor = rg),
    colGraphs = list(cell_cor = cg)
)
ge
#> class: GraphExperiment 
#> dim: 200 100 
#> metadata(0):
#> assays(1): counts
#> rownames(200): gene1 gene2 ... gene199 gene200
#> rowData names(3): pathway coding gene_cor__degree
#> colnames(100): cell1 cell2 ... cell99 cell100
#> colData names(1): cell_type
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
#> rowGraphs(1): gene_cor
#> colGraphs(1): cell_cor
```

If you’re familiar with `SummarizedExperiment` and
`SingleCellExperiment` objects, you will certainly recognize nearly
everything you see in `ge`. Compared to `SingleCellExperiment` objects,
the only difference here is in the last two rows, which indicate that
this object contains a `rowGraph` named ‘gene_cor’ and a `colGraph`
named ‘cell_cor’.

Importantly, since nodes of `rowGraphs`/`colGraphs` are always in sync
with `rownames`/`colnames`, **feature IDs in rownames and rowGraphs must
be the same**, and likewise for observation IDs in colnames and
colGraphs. For example, attempting to create a `GraphExperiment` object
with some features from `rownames` missing would lead to an error:

``` r

# Remove 'gene1' to 'gene10' from the rowGraph and try to recreate object
rg2 <- delete_vertices(rg, paste0("gene", 1:10))
GraphExperiment(
    assays = list(counts = mat),
    rowData = rdata,
    colData = cdata,
    rowGraphs = list(gene_cor = rg2)
)
#> Error in `validObject()`:
#> ! invalid class "GraphExperiment" object: 
#>     10 feature(s) in 'rownames' are missing from graph 'gene_cor'.
```

Alternatively, you can create a `GraphExperiment` object by coercing
from an existing `(Ranged)SummarizedExperiment` or
`SingleCellExperiment` object. For example:

``` r

# Coercing from `SummarizedExperiment`
se <- SummarizedExperiment(list(counts = mat))
ge1 <- as(se, "GraphExperiment")
ge1
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

Note that the `rowGraphs`/`colGraphs` containers are still there, but
empty.

To access the names of all graphs, you will use the
[`rowGraphNames()`](../reference/GraphExperiment-methods.md) and
[`colGraphNames()`](../reference/GraphExperiment-methods.md) functions.

``` r

# Get rowGraph names
rowGraphNames(ge)      # 'gene_cor'
#> [1] "gene_cor"
rowGraphNames(ge1)     # empty (NULL)
#> NULL

# Get colGraph names
colGraphNames(ge)      # 'cell_cor'
#> [1] "cell_cor"
colGraphNames(ge1)     # empty (NULL)
#> NULL
```

## Accessing `rowGraphs`/`colGraphs` and `rowData`/`colData` (a.k.a. ‘getters’)

To access graphs in `rowGraphs`/`colGraphs`, you can use one of the
following getter functions:

- `rowGraphs(x)`/`colGraphs(x)`: retrieves **all** (row/col)Graphs as a
  list of `igraph` objects.
- `rowGraph(x, i)`/`colGraph(x, i)`: retrieves only graph $`i`$ from the
  list. Note that $`i`$ can be a numeric scalar (index) or a character
  scalar (name).

The design here is equivalent to `assays()` versus `assay()` for
`SummarizedExperiment` objects.

``` r

# Get rowGraphs
rowGraphs(ge)
#> List of length 1
#> names(1): gene_cor

# Get colGraphs
colGraphs(ge)
#> List of length 1
#> names(1): cell_cor

# Get first rowGraph by index
rowGraph(ge, 1)
#> IGRAPH 8c54443 UNW- 200 20096 -- 
#> + attr: name (v/c), degree (v/n), pathway (v/c), coding (v/l), weight
#> | (e/n)
#> + edges from 8c54443 (vertex names):
#>  [1] gene1--gene1  gene1--gene2  gene1--gene3  gene1--gene4  gene1--gene5 
#>  [6] gene1--gene6  gene1--gene7  gene1--gene8  gene1--gene9  gene1--gene10
#> [11] gene1--gene11 gene1--gene12 gene1--gene13 gene1--gene14 gene1--gene15
#> [16] gene1--gene16 gene1--gene17 gene1--gene18 gene1--gene19 gene1--gene20
#> [21] gene1--gene21 gene1--gene22 gene1--gene23 gene1--gene24 gene1--gene25
#> [26] gene1--gene26 gene1--gene27 gene1--gene28 gene1--gene29 gene1--gene30
#> [31] gene1--gene31 gene1--gene32 gene1--gene33 gene1--gene34 gene1--gene35
#> + ... omitted several edges

# Get first rowGraph by index (alternative)
rowGraphs(ge)[[1]]
#> IGRAPH 8c54443 UNW- 200 20096 -- 
#> + attr: name (v/c), degree (v/n), pathway (v/c), coding (v/l), weight
#> | (e/n)
#> + edges from 8c54443 (vertex names):
#>  [1] gene1--gene1  gene1--gene2  gene1--gene3  gene1--gene4  gene1--gene5 
#>  [6] gene1--gene6  gene1--gene7  gene1--gene8  gene1--gene9  gene1--gene10
#> [11] gene1--gene11 gene1--gene12 gene1--gene13 gene1--gene14 gene1--gene15
#> [16] gene1--gene16 gene1--gene17 gene1--gene18 gene1--gene19 gene1--gene20
#> [21] gene1--gene21 gene1--gene22 gene1--gene23 gene1--gene24 gene1--gene25
#> [26] gene1--gene26 gene1--gene27 gene1--gene28 gene1--gene29 gene1--gene30
#> [31] gene1--gene31 gene1--gene32 gene1--gene33 gene1--gene34 gene1--gene35
#> + ... omitted several edges

# Get graph by name
rowGraph(ge, "gene_cor")
#> IGRAPH 8c54443 UNW- 200 20096 -- 
#> + attr: name (v/c), degree (v/n), pathway (v/c), coding (v/l), weight
#> | (e/n)
#> + edges from 8c54443 (vertex names):
#>  [1] gene1--gene1  gene1--gene2  gene1--gene3  gene1--gene4  gene1--gene5 
#>  [6] gene1--gene6  gene1--gene7  gene1--gene8  gene1--gene9  gene1--gene10
#> [11] gene1--gene11 gene1--gene12 gene1--gene13 gene1--gene14 gene1--gene15
#> [16] gene1--gene16 gene1--gene17 gene1--gene18 gene1--gene19 gene1--gene20
#> [21] gene1--gene21 gene1--gene22 gene1--gene23 gene1--gene24 gene1--gene25
#> [26] gene1--gene26 gene1--gene27 gene1--gene28 gene1--gene29 gene1--gene30
#> [31] gene1--gene31 gene1--gene32 gene1--gene33 gene1--gene34 gene1--gene35
#> + ... omitted several edges
```

Careful readers will notice that this `igraph` object has node
attributes that were not present in the original graph: ‘pathway’ and
‘coding’. This is because
[`rowGraphs()`](../reference/GraphExperiment-methods.md)/[`rowGraph()`](../reference/GraphExperiment-methods.md)
automatically extract `rowData` variables (if any) and add them to node
attributes.
[`colGraphs()`](../reference/GraphExperiment-methods.md)/[`colGraph()`](../reference/GraphExperiment-methods.md)
work in the same way (but with `colData`, of course). The same happens
in the other direction: the
[`rowData()`](../reference/GraphExperiment-methods.md)/[`colData()`](../reference/GraphExperiment-methods.md)
methods for `GraphExperiment` objects automatically add node attributes
(if any) to `rowData`/`colData`.

``` r

# `rowGraphs` and `rowData` are always in sync!
rowData(ge)
#> DataFrame with 200 rows and 3 columns
#>             pathway    coding gene_cor__degree
#>         <character> <logical>        <numeric>
#> gene1            P1     FALSE         3.125219
#> gene2            P1      TRUE         4.454433
#> gene3            P2      TRUE         0.873847
#> gene4            P1     FALSE         1.511270
#> gene5            P1      TRUE         0.568855
#> ...             ...       ...              ...
#> gene196          P1      TRUE         0.773387
#> gene197          P2     FALSE         0.857870
#> gene198          P1      TRUE         1.591484
#> gene199          P1      TRUE         2.324144
#> gene200          P1     FALSE        -0.887166

# `colGraphs` and `colData` too - yay!
colGraph(ge, 1) # note the `cell_type` attribute extracted from `colData`
#> IGRAPH 98a5462 UNW- 100 5050 -- 
#> + attr: name (v/c), cell_type (v/c), weight (e/n)
#> + edges from 98a5462 (vertex names):
#>  [1] cell1--cell1  cell1--cell2  cell1--cell3  cell1--cell4  cell1--cell5 
#>  [6] cell1--cell6  cell1--cell7  cell1--cell8  cell1--cell9  cell1--cell10
#> [11] cell1--cell11 cell1--cell12 cell1--cell13 cell1--cell14 cell1--cell15
#> [16] cell1--cell16 cell1--cell17 cell1--cell18 cell1--cell19 cell1--cell20
#> [21] cell1--cell21 cell1--cell22 cell1--cell23 cell1--cell24 cell1--cell25
#> [26] cell1--cell26 cell1--cell27 cell1--cell28 cell1--cell29 cell1--cell30
#> [31] cell1--cell31 cell1--cell32 cell1--cell33 cell1--cell34 cell1--cell35
#> [36] cell1--cell36 cell1--cell37 cell1--cell38 cell1--cell39 cell1--cell40
#> + ... omitted several edges
```

Variables ‘pathway’ and ‘coding’ were in the original data frame we used
as `rowData`, but variable ’gene_cor\_\_degree’ was added by extracting
the *degree* attribute of nodes in rowGraph `gene_cor`.

## Modifying `GraphExperiment` objects (a.k.a. ‘setters’)

Like in the `SummarizedExperiment` and `SingleCellExperiment` classes,
all getter methods specific to `GraphExperiment` objects have a
corresponding setter method. Such methods allow users to modify elements
by adding `<-` after the getter method. For example, to add or replace a
particular graph, you would use the `rowGraph<-`/`colGraph<-` method as
follows:

``` r

# Create a new rowGraph without correlations between -0.4 and 0.4
rg_filt <- rowGraph(ge, "gene_cor") |> 
    delete_vertex_attr("pathway") |>
    delete_vertex_attr("degree") |>
    delete_vertex_attr("coding")

todelete <- abs(E(rg_filt)$weight) <0.4
rg_filt <- delete_edges(rg_filt, which(todelete))
rg_filt
#> IGRAPH abfa81f UNW- 200 202 -- 
#> + attr: name (v/c), weight (e/n)
#> + edges from abfa81f (vertex names):
#>  [1] gene1 --gene1  gene2 --gene2  gene3 --gene3  gene4 --gene4  gene5 --gene5 
#>  [6] gene6 --gene6  gene7 --gene7  gene8 --gene8  gene9 --gene9  gene10--gene10
#> [11] gene11--gene11 gene12--gene12 gene13--gene13 gene14--gene14 gene15--gene15
#> [16] gene16--gene16 gene17--gene17 gene18--gene18 gene19--gene19 gene20--gene20
#> [21] gene21--gene21 gene22--gene22 gene23--gene23 gene24--gene24 gene25--gene25
#> [26] gene26--gene26 gene27--gene27 gene28--gene28 gene29--gene29 gene30--gene30
#> [31] gene31--gene31 gene32--gene32 gene33--gene33 gene34--gene34 gene35--gene35
#> [36] gene36--gene36 gene37--gene37 gene38--gene38 gene39--gene39 gene40--gene40
#> + ... omitted several edges

# Add filtered graph a new graph named `fcor`
rowGraph(ge, "filt_genecor") <- rg_filt
ge
#> class: GraphExperiment 
#> dim: 200 100 
#> metadata(0):
#> assays(1): counts
#> rownames(200): gene1 gene2 ... gene199 gene200
#> rowData names(3): pathway coding gene_cor__degree
#> colnames(100): cell1 cell2 ... cell99 cell100
#> colData names(1): cell_type
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
#> rowGraphs(2): gene_cor filt_genecor
#> colGraphs(1): cell_cor
```

If you’d like to replace all graphs at once, you could use the
`rowGraphs<-`/`colGraphs<-` setters. For example, let’s add a few graphs
to the `GraphExperiment` object we created before by coercing from
`SummarizedExperiment`:

``` r

# Taking a quick look (note: nothing in `rowGraphs`/`colGraphs`)
ge1
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

# Adding graphs from `ge`
rowGraphs(ge1) <- rowGraphs(ge)
colGraphs(ge1) <- colGraphs(ge)
ge1
#> class: GraphExperiment 
#> dim: 200 100 
#> metadata(0):
#> assays(1): counts
#> rownames(200): gene1 gene2 ... gene199 gene200
#> rowData names(5): gene_cor__degree gene_cor__pathway gene_cor__coding
#>   filt_genecor__pathway filt_genecor__coding
#> colnames(100): cell1 cell2 ... cell99 cell100
#> colData names(1): cell_cor__cell_type
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
#> rowGraphs(2): gene_cor filt_genecor
#> colGraphs(1): cell_cor
```

Lastly, you can also rename graphs by updating
`rowGraphNames`/`colGraphNames` as follows:

``` r

# Rename graphs
rowGraphNames(ge1) <- c("correlations", "correlations_filtered_0.4")
colGraphNames(ge1) <- c("cell_correlations")
ge1
#> class: GraphExperiment 
#> dim: 200 100 
#> metadata(0):
#> assays(1): counts
#> rownames(200): gene1 gene2 ... gene199 gene200
#> rowData names(5): correlations__degree correlations__pathway
#>   correlations__coding correlations_filtered_0.4__pathway
#>   correlations_filtered_0.4__coding
#> colnames(100): cell1 cell2 ... cell99 cell100
#> colData names(1): cell_correlations__cell_type
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
#> rowGraphs(2): correlations correlations_filtered_0.4
#> colGraphs(1): cell_correlations
```

## Subsetting `GraphExperiment` objects

In `SummarizedExperiment` objects, subsetting rows and columns (using
square brackets, `[`) automatically subsets `rowData` and `colData`
besides the assays. The same is true for `SingleCellExperiment` objects:
subsetting columns automatically subsets `colData` and `reducedDims`.

Since graphs in `GraphExperiment` objects are linked to rows and
columns, subsetting rows of a `GraphExperiment` object automatically
subsets rows of the `assays`, `rowData`, and all graphs in `rowGraphs`,
and subsetting columns automatically subsets columns of the `assays`,
`colData`, and graphs in `colGraphs`. For example:

``` r

# Subsetting `GraphExperiment` object
ge_subset <- ge[1:10, 1:10]

ge_subset
#> class: GraphExperiment 
#> dim: 10 10 
#> metadata(0):
#> assays(1): counts
#> rownames(10): gene1 gene2 ... gene9 gene10
#> rowData names(3): pathway coding gene_cor__degree
#> colnames(10): cell1 cell2 ... cell9 cell10
#> colData names(1): cell_type
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
#> rowGraphs(2): gene_cor filt_genecor
#> colGraphs(1): cell_cor
rowGraph(ge_subset, "gene_cor")
#> IGRAPH 2a24a7c UNW- 10 55 -- 
#> + attr: name (v/c), degree (v/n), pathway (v/c), coding (v/l), weight
#> | (e/n)
#> + edges from 2a24a7c (vertex names):
#>  [1] gene1--gene1 gene1--gene2 gene2--gene2 gene1--gene3 gene2--gene3
#>  [6] gene3--gene3 gene1--gene4 gene2--gene4 gene3--gene4 gene4--gene4
#> [11] gene1--gene5 gene2--gene5 gene3--gene5 gene4--gene5 gene5--gene5
#> [16] gene1--gene6 gene2--gene6 gene3--gene6 gene4--gene6 gene5--gene6
#> [21] gene6--gene6 gene1--gene7 gene2--gene7 gene3--gene7 gene4--gene7
#> [26] gene5--gene7 gene6--gene7 gene7--gene7 gene1--gene8 gene2--gene8
#> [31] gene3--gene8 gene4--gene8 gene5--gene8 gene6--gene8 gene7--gene8
#> + ... omitted several edges
```

## Session information

This document was created under the following conditions:

``` r

sessioninfo::session_info()
#> ─ Session info ───────────────────────────────────────────────────────────────
#>  setting  value
#>  version  R version 4.6.0 (2026-04-24)
#>  os       Ubuntu 24.04.4 LTS
#>  system   x86_64, linux-gnu
#>  ui       X11
#>  language en
#>  collate  en_US.UTF-8
#>  ctype    en_US.UTF-8
#>  tz       UTC
#>  date     2026-06-17
#>  pandoc   3.9.0.2 @ /usr/bin/ (via rmarkdown)
#>  quarto   1.9.38 @ /usr/local/bin/quarto
#> 
#> ─ Packages ───────────────────────────────────────────────────────────────────
#>  package              * version date (UTC) lib source
#>  abind                  1.4-8   2024-09-12 [1] RSPM (R 4.6.0)
#>  Biobase              * 2.73.1  2026-04-29 [1] Bioconductor 3.24 (R 4.6.0)
#>  BiocBaseUtils          1.15.1  2026-05-10 [1] Bioconductor 3.24 (R 4.6.0)
#>  BiocGenerics         * 0.59.7  2026-06-07 [1] Bioconductor 3.24 (R 4.6.0)
#>  BiocManager            1.30.27 2025-11-14 [1] RSPM (R 4.6.0)
#>  BiocStyle            * 2.41.0  2026-04-28 [1] Bioconductor 3.24 (R 4.6.0)
#>  bookdown               0.47    2026-06-16 [1] RSPM (R 4.6.0)
#>  bslib                  0.11.0  2026-05-16 [2] RSPM (R 4.6.0)
#>  cachem                 1.1.0   2024-05-16 [2] RSPM (R 4.6.0)
#>  cli                    3.6.6   2026-04-09 [2] RSPM (R 4.6.0)
#>  DelayedArray           0.39.3  2026-06-01 [1] Bioconductor 3.24 (R 4.6.0)
#>  desc                   1.4.3   2023-12-10 [2] RSPM (R 4.6.0)
#>  digest                 0.6.39  2025-11-19 [2] RSPM (R 4.6.0)
#>  evaluate               1.0.5   2025-08-27 [2] RSPM (R 4.6.0)
#>  fastmap                1.2.0   2024-05-15 [2] RSPM (R 4.6.0)
#>  fs                     2.1.0   2026-04-18 [2] RSPM (R 4.6.0)
#>  generics             * 0.1.4   2025-05-09 [1] RSPM (R 4.6.0)
#>  GenomicRanges        * 1.65.0  2026-04-28 [1] Bioconductor 3.24 (R 4.6.0)
#>  glue                   1.8.1   2026-04-17 [2] RSPM (R 4.6.0)
#>  GraphExperiment      * 1.1.2   2026-06-17 [1] Bioconductor
#>  htmltools              0.5.9   2025-12-04 [2] RSPM (R 4.6.0)
#>  htmlwidgets            1.6.4   2023-12-06 [2] RSPM (R 4.6.0)
#>  igraph               * 2.3.2   2026-05-29 [1] RSPM (R 4.6.0)
#>  IRanges              * 2.47.2  2026-06-01 [1] Bioconductor 3.24 (R 4.6.0)
#>  jquerylib              0.1.4   2021-04-26 [2] RSPM (R 4.6.0)
#>  jsonlite               2.0.0   2025-03-27 [2] RSPM (R 4.6.0)
#>  knitr                  1.51    2025-12-20 [2] RSPM (R 4.6.0)
#>  lattice                0.22-9  2026-02-09 [3] CRAN (R 4.6.0)
#>  lifecycle              1.0.5   2026-01-08 [2] RSPM (R 4.6.0)
#>  magrittr               2.0.5   2026-04-04 [2] RSPM (R 4.6.0)
#>  Matrix                 1.7-5   2026-03-21 [3] CRAN (R 4.6.0)
#>  MatrixGenerics       * 1.25.0  2026-04-28 [1] Bioconductor 3.24 (R 4.6.0)
#>  matrixStats          * 1.5.0   2025-01-07 [1] RSPM (R 4.6.0)
#>  otel                   0.2.0   2025-08-29 [2] RSPM (R 4.6.0)
#>  pillar                 1.11.1  2025-09-17 [2] RSPM (R 4.6.0)
#>  pkgconfig              2.0.3   2019-09-22 [2] RSPM (R 4.6.0)
#>  pkgdown                2.2.0   2025-11-06 [1] RSPM (R 4.6.0)
#>  R6                     2.6.1   2025-02-15 [2] RSPM (R 4.6.0)
#>  ragg                   1.5.2   2026-03-23 [2] RSPM (R 4.6.0)
#>  rlang                  1.2.0   2026-04-06 [2] RSPM (R 4.6.0)
#>  rmarkdown              2.31    2026-03-26 [1] RSPM (R 4.6.0)
#>  S4Arrays               1.13.0  2026-04-28 [1] Bioconductor 3.24 (R 4.6.0)
#>  S4Vectors            * 0.51.3  2026-06-01 [1] Bioconductor 3.24 (R 4.6.0)
#>  sass                   0.4.10  2025-04-11 [2] RSPM (R 4.6.0)
#>  Seqinfo              * 1.3.0   2026-04-28 [1] Bioconductor 3.24 (R 4.6.0)
#>  sessioninfo            1.2.4   2026-06-04 [2] RSPM (R 4.6.0)
#>  SingleCellExperiment * 1.35.1  2026-05-14 [1] Bioconductor 3.24 (R 4.6.0)
#>  SparseArray            1.13.2  2026-05-01 [1] Bioconductor 3.24 (R 4.6.0)
#>  SummarizedExperiment * 1.43.0  2026-04-28 [1] Bioconductor 3.24 (R 4.6.0)
#>  systemfonts            1.3.2   2026-03-05 [2] RSPM (R 4.6.0)
#>  textshaping            1.0.5   2026-03-06 [2] RSPM (R 4.6.0)
#>  vctrs                  0.7.3   2026-04-11 [2] RSPM (R 4.6.0)
#>  xfun                   0.58    2026-06-01 [2] RSPM (R 4.6.0)
#>  XVector                0.53.0  2026-04-28 [1] Bioconductor 3.24 (R 4.6.0)
#>  yaml                   2.3.12  2025-12-10 [2] RSPM (R 4.6.0)
#> 
#>  [1] /__w/_temp/Library
#>  [2] /usr/local/lib/R/site-library
#>  [3] /usr/local/lib/R/library
#>  * ── Packages attached to the search path.
#> 
#> ──────────────────────────────────────────────────────────────────────────────
```

## References

Amezquita, Robert, Aaron Lun, Etienne Becht, et al. 2020. “Orchestrating
Single-Cell Analysis with Bioconductor.” *Nature Methods* 17: 137–45.
<https://www.nature.com/articles/s41592-019-0654-x>.

[^1]: **Note on software design:** if you’re familiar with
    `SingleCellExperiment` objects, you probably know that it offers
    `rowPairs`/`colPairs` slots to store pairwise relationships between
    rows and columns of assays, respectively. In theory, some of the
    data stored in `rowGraphs`/`colGraphs` (of a `GraphExperiment`
    object) could be stored in `rowPairs`/`colPairs` (of a
    `SingleCellExperiment`). However, we chose to implement a dedicated
    slot with `igraph` objects to guarantee (i) seamless
    interoperability with other packages, given that `igraph` is the de
    facto standard class for graphs in R; and (ii) convenience in
    methods (e.g., subsetting, integration with `rowData`/`colData`,
    integration across multiple graphs, etc).

[^2]: **Tip:** in day-to-day single-cell RNA-seq analyses, researchers
    typically infer cell-cell graphs based on shared nearest-neighbors
    (SNN), which are then used to find clusters that can be mapped to
    cell types. Readers interested in this sort of graph can have a look
    at the `buildSNNGraph()` function from the
    *[scran](https://bioconductor.org/packages/3.24/scran)* package.
