# nice ggplot theme

nice theme with square border

## Usage

``` r
theme_rect()
```

## Examples

``` r
library(ggplot2)
p <- ggplot(mtcars, aes(x = wt, y = mpg)) +
  geom_point()
p + theme_rect()
```
