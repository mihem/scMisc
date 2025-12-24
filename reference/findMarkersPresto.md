# wrapper function for presto find markers

find significant DE genes using presto

## Usage

``` r
findMarkersPresto(
  ident1,
  ident2 = NULL,
  object,
  only_pos = FALSE,
  min_pct = 0.1,
  logfc_threshold = 0.25,
  assay = assay
)
```

## Arguments

- ident1:

  cell population 1

- ident2:

  cell population 2, if NULL use all (default: NULL)

- object:

  Seurat object

- only_pos:

  only return positive markers? (default: FALSE)

- min_pct:

  minimum fraction of cells in either two of the populations (default:
  0.1)

- logfc_threshold:

  minimum x-fold difference (default: 0.25)

- assay:

  which assay to use in DE testing (e.g. RNA or SCT)

## Value

data frame with significant DE genes arranged by log2FC

## Examples

``` r
library(Seurat)
findMarkersPresto(ident1 = "0", ident2 = "1", object = pbmc_small, assay = "RNA")
#> For a (much!) faster implementation of the Wilcoxon Rank Sum Test,
#> (default method for FindMarkers) please install the presto package
#> --------------------------------------------
#> install.packages('devtools')
#> devtools::install_github('immunogenomics/presto')
#> --------------------------------------------
#> After installation of presto, Seurat will automatically use the more 
#> efficient implementation (no further action necessary).
#> This message will be shown once per session
#>         gene avg_log2FC        p_val    p_val_adj pct.1 pct.2
#> 1        CD7  11.172163 2.291671e-05 5.270843e-03 0.528  0.00
#> 2        LCK   6.298467 1.179184e-04 2.712123e-02 0.500  0.04
#> 3       CCL5   5.198553 6.564657e-05 1.509871e-02 0.667  0.28
#> 4     NFKBIA  -1.223713 1.892762e-04 4.353352e-02 0.278  0.88
#> 5       PSAP  -1.375977 5.632116e-05 1.295387e-02 0.417  1.00
#> 6       SAT1  -1.473774 1.394154e-04 3.206554e-02 0.639  1.00
#> 7     FCER1G  -1.513784 6.370800e-06 1.465284e-03 0.389  1.00
#> 8     CARD16  -1.630200 8.720825e-05 2.005790e-02 0.139  0.72
#> 9     TYROBP  -1.714096 4.686298e-07 1.077848e-04 0.444  1.00
#> 10 LINC00936  -1.800795 4.973217e-05 1.143840e-02 0.083  0.60
#> 11     COTL1  -1.829941 2.124963e-05 4.887414e-03 0.500  0.96
#> 12   S100A11  -1.839553 1.235720e-06 2.842157e-04 0.333  1.00
#> 13    LGALS1  -1.857736 5.123692e-07 1.178449e-04 0.500  1.00
#> 14     GSTP1  -2.010925 1.100905e-05 2.532081e-03 0.333  0.88
#> 15      CTSS  -2.078919 8.095594e-07 1.861987e-04 0.361  0.96
#> 16      TSPO  -2.148317 8.902375e-06 2.047546e-03 0.333  0.88
#> 17      FCN1  -2.438790 5.418631e-07 1.246285e-04 0.167  0.88
#> 18      AIF1  -2.626718 8.106490e-09 1.864493e-06 0.194  1.00
#> 19    LGALS3  -2.639987 4.816986e-06 1.107907e-03 0.111  0.72
#> 20       CFP  -2.819252 6.044000e-06 1.390120e-03 0.083  0.68
#> 21    NUP214  -2.820961 5.481682e-05 1.260787e-02 0.028  0.48
#> 22       GRN  -2.829504 7.180762e-07 1.651575e-04 0.083  0.76
#> 23       BID  -2.875011 2.976826e-06 6.846699e-04 0.083  0.72
#> 24   HLA-DMA  -3.007053 2.689683e-05 6.186270e-03 0.056  0.56
#> 25      CTSB  -3.045623 1.993113e-05 4.584160e-03 0.111  0.68
#> 26     FCGRT  -3.078474 3.166989e-08 7.284074e-06 0.083  0.84
#> 27       HCK  -3.092237 6.751546e-06 1.552856e-03 0.028  0.56
#> 28      LST1  -3.099086 1.969061e-09 4.528840e-07 0.139  0.96
#> 29    LGALS2  -3.247194 9.650274e-06 2.219563e-03 0.056  0.60
#> 30    IFITM3  -3.332581 2.301085e-08 5.292495e-06 0.056  0.80
#> 31    RNF130  -3.338908 2.473344e-05 5.688692e-03 0.056  0.56
#> 32     IFI30  -3.387459 6.624184e-07 1.523562e-04 0.028  0.64
#> 33   HLA-DRA  -3.450292 1.268278e-08 2.917039e-06 0.417  0.92
#> 34  HLA-DQB1  -3.505836 5.620640e-06 1.292747e-03 0.028  0.56
#> 35  SERPINA1  -3.540765 1.886774e-07 4.339580e-05 0.056  0.72
#> 36       CFD  -3.659355 1.963076e-06 4.515076e-04 0.056  0.64
#> 37      CST3  -3.710886 4.469249e-11 1.027927e-08 0.306  1.00
#> 38      TYMP  -3.806179 1.702818e-11 3.916481e-09 0.111  1.00
#> 39    LRRC25  -3.923695 1.256935e-04 2.890952e-02 0.028  0.44
#> 40  HLA-DRB5  -4.122944 3.920170e-08 9.016391e-06 0.056  0.76
#> 41  HLA-DPA1  -4.197597 7.347202e-09 1.689856e-06 0.111  0.84
#> 42       LYZ  -4.476489 6.997602e-11 1.609449e-08 0.417  1.00
#> 43  HLA-DQA1  -4.711130 1.136284e-05 2.613452e-03 0.028  0.52
#> 44      IL1B  -4.725789 2.948315e-05 6.781125e-03 0.028  0.48
#> 45     MPEG1  -4.803084 1.834020e-04 4.218246e-02 0.028  0.40
#> 46      FPR1  -4.912697 2.460668e-05 5.659537e-03 0.028  0.48
#> 47   HLA-DMB  -4.920761 2.145323e-05 4.934243e-03 0.028  0.48
#> 48  HLA-DRB1  -4.999548 3.287672e-10 7.561646e-08 0.083  0.88
#> 49    S100A9  -5.125239 1.542308e-09 3.547309e-07 0.194  0.88
#> 50  HLA-DPB1  -5.325399 4.061018e-10 9.340340e-08 0.083  0.88
#> 51    S100A8  -6.076803 5.334985e-11 1.227047e-08 0.111  0.96
#> 52     SMCO4  -9.369438 1.739832e-06 4.001614e-04 0.000  0.52
#> 53     C5AR1  -9.594401 1.556284e-05 3.579453e-03 0.000  0.44
#> 54      CD14  -9.883809 1.739832e-06 4.001614e-04 0.000  0.52
#> 55    MS4A6A -10.301974 1.739832e-06 4.001614e-04 0.000  0.52
```
