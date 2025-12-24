# nice pheatmap

create and save a nice pheatmap in the folder `heatmap` using color
breaks

## Usage

``` r
pHeatmap(
  matrix,
  scale = "none",
  height = ceiling(nrow(matrix)/3),
  width = ceiling(ncol(matrix)/2),
  cellwidth = 10,
  cellheight = 10,
  treeheight_row = 10,
  treeheight_col = 10,
  fontsize = 10,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_row = NA,
  dir_output = "."
)
```

## Arguments

- matrix:

  numeric matrix of the values to be plotted

- scale:

  should the values be centered and scaled in row, column direction?
  Allowed: `row`, `column`, `none` (default: `none`)

- height:

  height of output plot (default: number of rows divided by three)

- width:

  width of output plot (default: number of columns divided by two)

- cellwidth:

  individual cell width (default: 10)

- cellheight:

  individual cell height (default: 10)

- treeheight_row:

  height of a tree for rows (default: 10)

- treeheight_col:

  height of a tree for columns (default: 10)

- fontsize:

  fontsize (default: 10)

- cluster_rows:

  cluster rows? (default: true)

- cluster_cols:

  cluster columns? (default: true)

- annotation_row:

  data frame that contains the annotations. Rows in the data and in the
  annotation are matched using row names. (default: NA)

- dir_output:

  directory to save the output plot (default: ".")

## Value

save heatmap to folder

## Examples

``` r
matrix <- matrix(rnorm(100), nrow = 10, ncol = 10)
rownames(matrix) <- paste0("Gene", 1:10)
colnames(matrix) <- paste0("Sample", 1:10)
pHeatmap(matrix, scale = "row", dir_output = ".")

#> agg_record_1140350788 
#>                     2 
unlink("hm_matrix.pdf")
```
