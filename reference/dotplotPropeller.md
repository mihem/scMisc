# plot propeller results in a barplot

The function creates a barplot of the propeller results and saves the
plot

## Usage

``` r
dotplotPropeller(
  data,
  color,
  filename,
  width = 5,
  height = 5,
  dir_output = "."
)
```

## Arguments

- data:

  A dataframe containing the results from propeller calculation

- color:

  A vector of colors for the clusters in the plot

- filename:

  A character representing the file name of the plot

- width:

  The width of the plot

- height:

  The height of the plot

- dir_output:

  directory to save the output plot (default: ".")

## Value

create and save propeller abundance barplot

## Examples

``` r
propeller_data <- data.frame(
     cluster = c("Cluster1", "Cluster2", "Cluster3"),
    log2ratio = c(1.5, -2.0, 0.5)
)
color <- c("Cluster1" = "blue", "Cluster2" = "red", "Cluster3" = "green")
dotplotPropeller(data = propeller_data,
 color = color,
 filename = "test_propeller_dotplot",
 width = 5,
 height = 5,
 dir_output = ".")

unlink("propeller_dotplot_test_propeller_dotplot.pdf")
```
