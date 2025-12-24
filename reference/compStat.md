# significance for boxplot

calculcate significance for ggsignif

## Usage

``` r
compStat(x_var, group, data, paired)
```

## Arguments

- x_var:

  numeric variable (x in the formula `x ~ group`)

- group:

  grouping variable (group in the formula `x ~ group`)

- data:

  data frame containing the variables

- paired:

  logical value. Do you want to do perform pair testing?

## Value

data frame with adjusted signficance values

## Examples

``` r
set.seed(123)
data <- data.frame(
  sample = c(paste0("CSF_P0", 1:9)),
  type = c(rep("control", 4), rep("AIE", 5)),
  cluster1 = c(runif(4, 0, 1), runif(5, 99, 100)),
  cluster2 = c(runif(4, 0, 70), runif(5, 20, 100))
)
compStat(x_var = c("cluster1", "cluster2"), group = "type", data = data, paired = FALSE)
#> # A tibble: 1 × 6
#>   .y.      group1 group2       p  p.adj p.adj.signif
#>   <chr>    <chr>  <chr>    <dbl>  <dbl> <chr>       
#> 1 cluster1 AIE    control 0.0159 0.0318 *           
```
