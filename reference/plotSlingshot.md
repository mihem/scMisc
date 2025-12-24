# plot slingshot results

The function creates a colored UMAP plot and curves split by lineage

## Usage

``` r
plotSlingshot(object, lineage, pt, curves)
```

## Arguments

- object:

  A Seurat object

- lineage:

  A character string indicating the lineage of interest.

- pt:

  A numeric vector of pseudotime values.

- curves:

  A dataframe containing the curves.

## Value

A ggplot object.

## Examples

``` r
library(Seurat)
set.seed(123)
pbmc_small$lineage <- sample(c("Lineage1", "Lineage2"), ncol(pbmc_small), replace = TRUE)
pbmc_small$umap <-
CreateDimReducObject(
 embeddings =
   Embeddings(
     pbmc_small,
     reduction = "tsne"
   ),
 key = "UMAP_",
 assay = "RNA"
)
curves <- data.frame(
    UMAP_1 = runif(ncol(pbmc_small), min = -10, max = 10),
    UMAP_2 = runif(ncol(pbmc_small), min = -10, max = 10),
    Lineage = sample(c("Lineage1", "Lineage2"), ncol(pbmc_small), replace = TRUE)
)
pt <- matrix(runif(ncol(pbmc_small) * 2), ncol = 2)
colnames(pt) <- c("Lineage1", "Lineage2")
plot <- plotSlingshot(object = pbmc_small, lineage = "Lineage1", pt = pt, curves = curves)
```
