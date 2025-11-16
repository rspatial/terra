# Histogram

Create a histogram of the values of a SpatRaster. For large datasets a
sample of `maxcell` is used.

## Usage

``` r
# S4 method for class 'SpatRaster'
hist(x, layer, maxcell=1000000, plot=TRUE, maxnl=16, main, ...)
```

## Arguments

- x:

  SpatRaster

- layer:

  positive integer or character to indicate layer numbers (or names). If
  missing, all layers up to `maxnl` are used

- maxcell:

  integer. To regularly sample very large objects

- plot:

  logical. Plot the histogram or only return the histogram values

- maxnl:

  positive integer. The maximum number of layers to use. Ignored if
  `layer` is not missing

- main:

  character. Main title(s) for the plot. Default is the value of
  [`names`](https://rspatial.github.io/terra/reference/names.md)

- ...:

  additional arguments. See
  [`hist`](https://rdrr.io/r/graphics/hist.html)

## Value

This function is principally used for plotting a histogram, but it also
returns an object of class "histogram" (invisibly if `plot=TRUE`).

## See also

[`pairs`](https://rspatial.github.io/terra/reference/pairs.md)`, `[`boxplot`](https://rspatial.github.io/terra/reference/boxplot.md)

## Examples

``` r
r1 <- r2 <- rast(nrows=50, ncols=50)
values(r1) <- runif(ncell(r1))
values(r2) <- runif(ncell(r1))
rs <- r1 + r2
rp <- r1 * r2

opar <- par(no.readonly =TRUE)
par(mfrow=c(2,2))
plot(rs, main='sum')
plot(rp, main='product')
hist(rs)
a <- hist(rp)

a
#> $breaks
#>  [1] 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0
#> 
#> $counts
#>  [1] 815 485 350 251 188 160 117  75  42  17
#> 
#> $density
#>  [1] 3.260 1.940 1.400 1.004 0.752 0.640 0.468 0.300 0.168 0.068
#> 
#> $mids
#>  [1] 0.05 0.15 0.25 0.35 0.45 0.55 0.65 0.75 0.85 0.95
#> 
#> $xname
#> [1] "lyr.1"
#> 
#> $equidist
#> [1] TRUE
#> 
#> attr(,"class")
#> [1] "histogram"
x <- c(rs, rp, sqrt(rs))
hist(x)
par(opar)
```
