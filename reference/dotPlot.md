# nice Seurat dot plot

create and save a nice Seurat dot plot

## Usage

``` r
dotPlot(
  path,
  object,
  par,
  dot_min,
  scale = TRUE,
  ortho = "none",
  width = 10,
  height = 10,
  dir_output = "."
)
```

## Arguments

- path:

  path to markers.csv

- object:

  Seurat object

- par:

  column name in markers.csv

- dot_min:

  minimal dot size

- scale:

  should the values be scaled? (default: TRUE)

- ortho:

  convert to orthologues? Allowed values: `none`, `mouse2human` or
  `human2mouse`

- width:

  width of output plot (default: 10)

- height:

  height of output plot (default: 10)

- dir_output:

  directory to save the output plot (default: ".")

## Value

save dot plot

## Examples

``` r
library(Seurat)
markers <- data.frame(B = c("MS4A1", "CD79A"))
write.csv(markers, "markers.csv")
dotPlot(
  path = "markers.csv",
  object = pbmc_small,
  par = "B",
  dot_min = 0.1,
  dir_output = "."
)
#> New names:
#> • `` -> `...1`
#> Rows: 2 Columns: 2
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: ","
#> chr (1): B
#> dbl (1): ...1
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
#> no genes were converted
#> Warning: Scaling data with a low number of groups may produce misleading results
#> Scale for colour is already present.
#> Adding another scale for colour, which will replace the existing scale.
#> Scale for size is already present.
#> Adding another scale for size, which will replace the existing scale.
#> Warning: Removed 3 rows containing missing values or values outside the scale range
#> (`geom_point()`).
#> Warning: Removed 3 rows containing missing values or values outside the scale range
#> (`geom_point()`).

unlink("markers.csv")
unlink("dp_pbmc_small_B.pdf")
```
