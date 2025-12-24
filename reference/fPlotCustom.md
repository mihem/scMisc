# Seurat feature plot

create and save a Seurat feature plot

## Usage

``` r
fPlotCustom(
  object,
  markers,
  par,
  reduction,
  width = 16,
  height = ceiling(length(genes_found)/4) * 3,
  dir_output = "."
)
```

## Arguments

- object:

  Seurat object

- markers:

  a data frame with a column called `cell_source` that represents the
  cell population and its source and a column `gene`

- par:

  a character string representing the cell_source to plot

- reduction:

  a character string specifying the dimension reduction

- width:

  width of output plot (default: 16)

- height:

  height of output plot (default: length of genes divided by four,
  ceiling, times three)

- dir_output:

  directory to save the output plot (default: ".")

## Value

save feature plot

## Examples

``` r
library(Seurat)
markers <- data.frame(cell_source = c("B", "B"), gene = c("MS4A1", "CD79A"))
fPlotCustom(
  object = pbmc_small,
  markers = markers,
  par = "B",
  reduction = "tsne"
)

unlink("fp_pbmc_small_B.png")
```
