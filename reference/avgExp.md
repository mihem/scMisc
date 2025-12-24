# wrapper function for Seurats AverageExpression

create average expression matrix with Seurat based on `markers.csv`

## Usage

``` r
avgExp(par, object, assay, slot, ortho = "none")
```

## Arguments

- par:

  column name in markers.csv

- object:

  Seurat object

- assay:

  which assay to use

- slot:

  which slot to use

- ortho:

  convert to orthologues? Allowed values: `none`, `mouse2human` or
  `human2mouse`

## Value

matrix with genes as rows and identity clases as columns

## Examples

``` r
library(Seurat)
markers <- data.frame(Bc = c("CD19", "MS4A1", "CD79B"))
write.csv(markers, "markers.csv")
avgExp("Bc", object = pbmc_small, assay = "RNA", slot = "data")
#> New names:
#> • `` -> `...1`
#> Rows: 3 Columns: 2
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: ","
#> chr (1): Bc
#> dbl (1): ...1
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
#> As of Seurat v5, we recommend using AggregateExpression to perform pseudo-bulk analysis.
#> This message is displayed once per session.
#> First group.by variable `ident` starts with a number, appending `g` to ensure valid variable names
#> This message is displayed once every 8 hours.
#> 3 x 3 sparse Matrix of class "dgCMatrix"
#>             g0        g1        g2
#> CD19   .        .         33.19047
#> MS4A1  .        2.083443 171.61521
#> CD79B 10.81466 17.548842 152.13441
unlink("markers.csv")
```
