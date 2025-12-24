# clustifyr UMAP heatmap wrap

create and save clustifyr UMAP and heatmap

## Usage

``` r
clustifyFun(
  my_object,
  my_annotations,
  ref_object,
  ref_annotations,
  format,
  ortho = "none",
  filter = FALSE,
  regex = "",
  umap = FALSE,
  width = 10,
  height = 7,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  cellwidth = 10,
  cellheight = 10,
  fontsize = 10,
  color,
  name
)
```

## Arguments

- my_object:

  query input object

- my_annotations:

  annotations of query input object

- ref_object:

  reference object

- ref_annotations:

  reference annotations

- format:

  format of reference object: possible values are seurat, sc or matrix

- ortho:

  convert to orthologues? Allowed values: `none`, `mouse2human` or
  `human2mouse` (default: `none`)

- filter:

  should the reference be filtered? Allowed values: TRUE, FALSE
  (default: FALSE)

- regex:

  if filter: specify filter term

- umap:

  also create UMAP? Allowed values: TRUE, FALSE (default FALSE)

- width:

  width of heatmap plot (default: 10)

- height:

  height of heatmap plot (default: 7)

- cluster_rows:

  cluster rows? (default: true)

- cluster_cols:

  cluster columns? (default: true)

- cellwidth:

  individual cell width (default: 10)

- cellheight:

  individual cell height (default: 10)

- fontsize:

  fontsize (default: 10)

- color:

  color palette

- name:

  fontsize (default: 10)

## Value

save umap/heatmap plot to folder `/results/clustifyr/`

## Examples

``` r
if (FALSE) { # \dontrun{
clustifyFun(
    my_object = tcells,
    my_annotations = "cluster",
    ref_object = projectil_til,
    ref_annotations = "functional.cluster",
    format = "seurat",
    ortho = "none",
    height = 25,
    width = 10,
    umap = TRUE,
    color = col_vector,
    name = "projectil_til"
)
} # }
```
