# GraphExperiment

**GraphExperiment** provides users and developers with infrastructure to
store quantitative data from omics assays along with networks
representing how assay features (e.g., genes, proteins, metabolites)
and/or observations (e.g., samples, cells, species) interact with each
other. The `GraphExperiment` S4 class extends `SingleCellExperiment` to
include support for `igraph` objects associated with features (stored in
`rowGraphs`) and observations (stored in `colGraphs`).

The figure belows summarizes the `GraphExperiment` class in comparison
to the `SingleCellExperiment` class.

![](articles/GraphExperiment.png)

## Installation instructions

Get the latest stable `R` release from
[CRAN](http://cran.r-project.org/). Then install `GraphExperiment` from
[Bioconductor](http://bioconductor.org/) using the following code:

``` r

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("GraphExperiment")
```

## Citation

Below is the citation output from using `citation('GraphExperiment')` in
R. Please run this yourself to check for any updates on how to cite
**GraphExperiment**.

``` r

print(citation('GraphExperiment'), bibtex = TRUE)
#> To cite package 'GraphExperiment' in publications use:
#> 
#>   Almeida-Silva F (2026). _GraphExperiment: S4 Class for Quantitative
#>   Data and Associated Networks_. R package version 1.1.2,
#>   <https://github.com/almeidasilvaf/GraphExperiment>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {GraphExperiment: S4 Class for Quantitative Data and Associated Networks},
#>     author = {Fabricio Almeida-Silva},
#>     year = {2026},
#>     note = {R package version 1.1.2},
#>     url = {https://github.com/almeidasilvaf/GraphExperiment},
#>   }
```

## Code of Conduct

Please note that the `GraphExperiment` project is released with a
[Contributor Code of
Conduct](http://bioconductor.org/about/code-of-conduct/). By
contributing to this project, you agree to abide by its terms.
