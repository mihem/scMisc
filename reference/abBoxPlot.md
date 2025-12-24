# abundance boxplot plot

create and save an abundance boxplot in the folder `abundance`

## Usage

``` r
abBoxPlot(
  object,
  cluster_idents,
  sample,
  cluster_order,
  group_by,
  group_order,
  color,
  width = 9,
  height = ceiling(length(unique(object@meta.data[[cluster_idents]]))/4) * 3,
  paired = FALSE,
  number_of_tests,
  dir_output = "."
)
```

## Arguments

- object:

  Seurat object

- cluster_idents:

  variable in meta data with cluster names

- sample:

  variable in meta data for each sample

- cluster_order:

  vector determining the order of the clusters

- group_by:

  variable in meta data that categorize samples in groups

- group_order:

  vector determining the order of the samples

- color:

  color palette

- width:

  width of output plot (default: 9)

- height:

  height of output plot (default: length of cluster_idents divided by
  four, ceiling, times three)

- paired:

  logical indicating whether you want a paired test (default FALSE)

- number_of_tests:

  number of tests to be performed

- dir_output:

  directory to save the output plot (default: ".")

## Value

save abundance box plot

## Examples

``` r
library(Seurat)
set.seed(123)
pbmc_small$cluster <-
  sample(c("Cluster1", "Cluster2"),
    ncol(pbmc_small),
    replace = TRUE
  )
pbmc_small$sample <-
  sample(c("CSF_P01", "CSF_P02", "CSF_P03", "CSF_P04"),
    ncol(pbmc_small),
    replace = TRUE
  )
lookup <-
  data.frame(
    sample = c("CSF_P01", "CSF_P02", "CSF_P03", "CSF_P04"),
    AIE_type = c(rep("control", 2), rep("CASPR2", 2))
  )
pbmc_small@meta.data <-
  pbmc_small@meta.data |>
  tibble::rownames_to_column("barcode") |>
  dplyr::left_join(lookup, by = "sample") |>
  tibble::column_to_rownames("barcode")
abBoxPlot(
  object = pbmc_small,
  cluster_idents = "cluster",
  sample = "sample",
  cluster_order = c("Cluster1", "Cluster2"),
  group_by = "AIE_type",
  group_order = c("control", "CASPR2"),
  color = c("control" = "blue", "CASPR2" = "red"),
  width = 9,
  height = 6,
  paired = FALSE,
  number_of_tests = 3,
  dir_output = "."
)
#> Joining with `by = join_by(sample)`
#> Warning: Computation failed in `stat_signif()`.
#> Caused by error in `$<-.data.frame`:
#> ! replacement has 1 row, data has 0
#> `geom_line()`: Each group consists of only one observation.
#> ℹ Do you need to adjust the group aesthetic?
#> Warning: Computation failed in `stat_signif()`.
#> Caused by error in `$<-.data.frame`:
#> ! replacement has 1 row, data has 0
#> `geom_line()`: Each group consists of only one observation.
#> ℹ Do you need to adjust the group aesthetic?
#> Warning: Computation failed in `stat_signif()`.
#> Caused by error in `$<-.data.frame`:
#> ! replacement has 1 row, data has 0
#> `geom_line()`: Each group consists of only one observation.
#> ℹ Do you need to adjust the group aesthetic?
#> Warning: Computation failed in `stat_signif()`.
#> Caused by error in `$<-.data.frame`:
#> ! replacement has 1 row, data has 0
#> `geom_line()`: Each group consists of only one observation.
#> ℹ Do you need to adjust the group aesthetic?

unlink("boxplot_cluster_pbmc_small_AIE_type.pdf")
```
