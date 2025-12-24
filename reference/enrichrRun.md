# enrichment analysis

wrapper function to perform enrichment analysis with enrichr and save
results

## Usage

``` r
enrichrRun(
  sheet,
  dir_input = ".",
  filename,
  dir_output = ".",
  dbs,
  fc_thresh = 1,
  p_thresh = 0.001,
  remove_rp_mt
)
```

## Arguments

- sheet:

  sheet name in excel file

- dir_input:

  path to the directory where the deg file is located (default: .)

- filename:

  name of deg file file (should be .xlsx) without extension

- dir_output:

  path to the directory where the results will be saved (default: .)

- dbs:

  name of the enrichr libraries

- fc_thresh:

  log fc threshold (default 1)

- p_thresh:

  p value threshold (default 0.001)

- remove_rp_mt:

  remove ribosomal and mitochondrial genes? (boolean value)

## Value

save enrichrment analysis in `dir_output` in two excel sheets (one for
positive and one for negative)

## Examples

``` r
library(Seurat)
library(enrichR)
#> Welcome to enrichR
#> Checking connections ... 
#> Enrichr ... 
#> Connection is Live!
#> FlyEnrichr ... 
#> Connection is Live!
#> WormEnrichr ... 
#> Connection is Live!
#> YeastEnrichr ... 
#> Connection is Live!
#> FishEnrichr ... 
#> Connection is Live!
#> OxEnrichr ... 
#> Connection is Live!
set.seed(123)
deg_data <- data.frame(
  gene = rownames(pbmc_small),
  avg_log2FC = rnorm(length(rownames(pbmc_small)), mean = 0, sd = 3),
  p_val_adj = runif(length(rownames(pbmc_small)), min = 0, max = 0.01)
)

writexl::write_xlsx(deg_data, path = "test_de.xlsx")
enrichrRun(
  sheet = "Sheet1",
  dir_input = ".",
  filename = "test_de",
  dbs = c("KEGG_2019_Human"),
  fc_thresh = 1,
  p_thresh = 0.001,
  remove_rp_mt = FALSE
)
#> Uploading data to Enrichr... Done.
#>   Querying KEGG_2019_Human... Done.
#> Parsing results... Done.
#> Uploading data to Enrichr... Done.
#>   Querying KEGG_2019_Human... Done.
#> Parsing results... Done.
unlink("test_de.xlsx")
unlink("enrichr_test_de_neg_Sheet1.xlsx")
unlink("enrichr_test_de_pos_Sheet1.xlsx")
```
