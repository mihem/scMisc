# plot PCA of cluster abundances

The function creates a PCA plot of the cluster abundances and saves the
plot in results/abundance folder

## Usage

``` r
pcaSeurat(
  object,
  cluster,
  sample,
  condition,
  width = 20,
  height = 5,
  dir_output = "."
)
```

## Arguments

- object:

  A Seurat object

- cluster:

  A character string indicating the cluster column in the metadata of
  the Seurat object

- sample:

  A character string indicating the sample column in the metadata of the
  Seurat object

- condition:

  A character string indicating the condition column in the metadata of
  the Seurat object

- width:

  The width of the plot

- height:

  The height of the plot

- dir_output:

  directory to save the output plot (default: ".")

## Value

Create and save PCA plot of Seurat cluster abundances.

## Examples

``` r
library(Seurat)
# Setup example data
set.seed(123)
pbmc_small$cluster <- sample(c("Cluster1", "Cluster2"), ncol(pbmc_small), replace = TRUE)
pbmc_small$sample <- sample(c("CSF_P01", "CSF_P02", "CSF_P03", "CSF_P04"), 
                           ncol(pbmc_small), replace = TRUE)

# Create lookup table for conditions
lookup <- data.frame(
  sample = c("CSF_P01", "CSF_P02", "CSF_P03", "CSF_P04"),
  condition = c(rep("control", 2), rep("treatment", 2))
)

# Add condition information to metadata
pbmc_small@meta.data <- pbmc_small@meta.data |>
  tibble::rownames_to_column("barcode") |>
  dplyr::left_join(lookup, by = "sample") |>
  tibble::column_to_rownames("barcode")

# Generate PCA plot
pcaSeurat(
  object = pbmc_small,
  cluster = "cluster",
  sample = "sample", 
  condition = "condition",
  width = 20,
  height = 5
)
#> Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
#> ℹ Please use `linewidth` instead.
#> ℹ The deprecated feature was likely used in the ggpubr package.
#>   Please report the issue at <https://github.com/kassambara/ggpubr/issues>.
#> Warning: `aes_string()` was deprecated in ggplot2 3.0.0.
#> ℹ Please use tidy evaluation idioms with `aes()`.
#> ℹ See also `vignette("ggplot2-in-packages")` for more information.
#> ℹ The deprecated feature was likely used in the factoextra package.
#>   Please report the issue at <https://github.com/kassambara/factoextra/issues>.
#> Joining with `by = join_by(cluster)`
#> Warning: Ignoring empty aesthetic: `width`.
#> Ignoring unknown labels:
#> • linetype : "group"
#> • shape : "group"
#> Warning: Ignoring empty aesthetic: `width`.
#> Ignoring unknown labels:
#> • linetype : "group"
#> • shape : "group"

unlink("pbmc_small_condition_cluster.pdf")
```
