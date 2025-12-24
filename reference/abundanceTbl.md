# wrapper function to calculate absolute and relative abundance

wrapper function to calculate absolute and relative abundance

## Usage

``` r
abundanceTbl(object, row_var, col_var, dir_output = ".")
```

## Arguments

- object:

  Seurat object

- row_var:

  variable in meta data that will represent the rows

- col_var:

  variable in meta data that will represent the columns

- dir_output:

  target directory to save the results (default: .)

## Examples

``` r
library(Seurat)
abundanceTbl(pbmc_small, row_var = "groups", col_var = "letter.idents")
unlink("abundance_tbl_pbmc_small_letter.idents.xlsx")
```
