# Load CellBender h5 matrices

Extract sparse matrix with corrected counts from CellBender h5 output
file.

## Usage

``` r
ReadCellBender_h5(file_name, use.names = TRUE, unique.features = TRUE)
```

## Arguments

- file_name:

  Path to h5 file.

- use.names:

  Label row names with feature names rather than ID numbers (default
  TRUE).

- unique.features:

  Make feature names unique (default TRUE).

## Value

sparse matrix

## References

Code used in function has been modified from
[`Seurat::Read10X_h5`](https://satijalab.org/seurat/reference/Read10X_h5.html)
function of Seurat package <https://github.com/satijalab/seurat>
(License: GPL-3).

## Examples

``` r
file_name <- system.file("extdata", "raw_feature_bc_matrix_filtered.h5", package = "scMisc")
mat <- ReadCellBender_h5(file_name)
```
