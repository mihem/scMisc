# abundance stacked bar plot

create and save an abundance stacked bar plot in the folder `abundance`

## Usage

``` r
stackedPlot(
  object,
  x_axis,
  y_axis,
  x_order,
  y_order,
  color,
  width,
  height = 10,
  dir_output = "."
)
```

## Arguments

- object:

  Seurat object

- x_axis:

  variable in meta data that is used for the y axis

- y_axis:

  variable in meta data that is used for the x axis

- x_order:

  vector determining the order of the x axis

- y_order:

  vector determining the order of the y axis

- color:

  color palette

- width:

  width of output plot (default: 10)

- height:

  height of output plot

- dir_output:

  directory to save the output plot (default: ".")

## Value

save stacked abundance barplot

## Examples

``` r
if (FALSE) { # \dontrun{
stackedPlot(
  object = sc_merge,
  x_axis = "pool",
  y_axis = "cluster",
  x_order = unique(sc_merge$pool),
  y_order = cluster_order,
  color = col_vector,
  width = 4
)
} # }
```
