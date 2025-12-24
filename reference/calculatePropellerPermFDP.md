# Calculate permFDP-adjusted threshold for propeller results

Internal helper function to calculate permFDP-adjusted significance
threshold for propeller differential abundance testing

## Usage

``` r
calculatePropellerPermFDP(
  prop_data,
  seu_obj,
  sample_col,
  cluster_col,
  meta_col,
  cond1,
  cond2,
  fdr_threshold,
  n_perms = 1000
)
```

## Arguments

- prop_data:

  A dataframe containing propeller results with P.Value column

- seu_obj:

  A Seurat object subset to the two conditions being compared

- sample_col:

  A character representing the name of the sample column

- cluster_col:

  A character representing the name of the cluster column

- meta_col:

  A character representing the name of the condition column

- cond1:

  A character representing the first condition

- cond2:

  A character representing the second condition

- fdr_threshold:

  The FDR threshold to use for permFDP adjustment

- n_perms:

  Number of permutations for permFDP adjustment (default: 1000)

## Value

A dataframe with added permFDP_sig and permFDP_threshold columns
