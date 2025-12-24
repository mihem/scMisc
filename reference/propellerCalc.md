# propellerCalc

perform a propeller t-test on a Seurat object

## Usage

``` r
propellerCalc(
  seu_obj1,
  condition1,
  condition2,
  cluster_col,
  meta_col,
  lookup,
  sample_col,
  formula,
  min_cells = 30,
  adjustment_method = "BH",
  fdr_threshold = NULL,
  n_perms = 1000
)
```

## Arguments

- seu_obj1:

  A seurat object

- condition1:

  A character representing the first condition

- condition2:

  A character representing the second condition

- cluster_col:

  A character representing the name of the column which contains cluster
  information

- meta_col:

  A character representing the name of the column which contains meta
  information, should be in Seurat object and in lookup table

- lookup:

  A dataframe that contains sample information

- sample_col:

  A character representing the name of the column in seu_obj1 which
  contains sample information

- formula:

  A linear model that should be used for the design matrix

- min_cells:

  A numeric value indicating the minimum number of cells in a cluster
  that should be included in the analysis

- adjustment_method:

  A character indicating the FDR adjustment method: "BH" (default) or
  "permFDP"

- fdr_threshold:

  The FDR threshold to use for permFDP adjustment (default: NULL)

- n_perms:

  Number of permutations for permFDP adjustment (default: 1000)

## Value

A dataframe containing the results of the t-test

## Examples

``` r
set.seed(123)
library(Seurat)
pbmc_small$condition <- factor(sample(c("diseaseA", "diseaseB"), nrow(pbmc_small), replace = TRUE))
pbmc_small$cluster <- Idents(pbmc_small)
pbmc_small$patient <- rep(paste0("P", 0:9), each = 8, length.out = ncol(pbmc_small))
lookup <- data.frame(
  patient = paste0("P", 0:9),
  condition = sample(c("diseaseA", "diseaseB"), 10, replace = TRUE)
)
# Standard BH adjustment
propellerCalc(
  seu_obj1 = pbmc_small,
  condition1 = "diseaseA",
  condition2 = "diseaseB",
  cluster_col = "cluster",
  meta_col = "condition",
  lookup = lookup,
  sample_col = "patient",
  formula = "~ 0 + condition",
  min_cells = 30
)
#> Performing logit transformation of proportions
#> # A tibble: 1 × 9
#>   cluster PropMean.conditiondiseaseA PropMean.conditiondi…¹ PropRatio Tstatistic
#>   <chr>                        <dbl>                  <dbl>     <dbl>      <dbl>
#> 1 0                            0.667                  0.357      1.87       1.28
#> # ℹ abbreviated name: ¹​PropMean.conditiondiseaseB
#> # ℹ 4 more variables: P.Value <dbl>, FDR <dbl>, log2ratio <dbl>, FDR_log <dbl>

# With permFDP adjustment (use lower min_cells to have multiple clusters)
propellerCalc(
  seu_obj1 = pbmc_small,
  condition1 = "diseaseA",
  condition2 = "diseaseB",
  cluster_col = "cluster",
  meta_col = "condition",
  lookup = lookup,
  sample_col = "patient",
  formula = "~ 0 + condition",
  min_cells = 5,
  adjustment_method = "permFDP",
  fdr_threshold = 0.1,
  n_perms = 100
)
#> Performing logit transformation of proportions
#> Running permFDP with 100 permutations for diseaseA vs diseaseB (3 clusters, 10 samples)...
#>   -> Corrected threshold: 0.0744
#> # A tibble: 3 × 11
#>   cluster PropMean.conditiondiseaseA PropMean.conditiondi…¹ PropRatio Tstatistic
#>   <chr>                        <dbl>                  <dbl>     <dbl>      <dbl>
#> 1 2                           0.0417                  0.321     0.130    -1.49  
#> 2 0                           0.667                   0.357     1.87      1.28  
#> 3 1                           0.292                   0.321     0.907    -0.0823
#> # ℹ abbreviated name: ¹​PropMean.conditiondiseaseB
#> # ℹ 6 more variables: P.Value <dbl>, FDR <dbl>, log2ratio <dbl>, FDR_log <dbl>,
#> #   permFDP_sig <lgl>, permFDP_threshold <dbl>
```
