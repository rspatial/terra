# summary

Compute summary statistics (min, max, mean, and quartiles) for
SpatRaster using base [`summary`](https://rdrr.io/r/base/summary.html)
method. A sample is used for very large files.

For single or other statistics see
[`Summary-methods`](https://rspatial.github.io/terra/reference/summarize-generics.md),
[`global`](https://rspatial.github.io/terra/reference/global.md), and
[`quantile`](https://rspatial.github.io/terra/reference/quantile.md)

## Usage

``` r
# S4 method for class 'SpatRaster'
summary(object, size=100000, warn=TRUE, ...)

# S4 method for class 'SpatVector'
summary(object, ...)
```

## Arguments

- object:

  SpatRaster or SpatVector

- size:

  positive integer. Size of a regular sample used for large datasets
  (see
  [`spatSample`](https://rspatial.github.io/terra/reference/sample.md))

- warn:

  logical. If `TRUE` a warning is given if a sample is used

- ...:

  additional arguments passed on to the base
  [`summary`](https://rdrr.io/r/base/summary.html)` method`

## Value

matrix with (an estimate of) the median, minimum and maximum values, the
first and third quartiles, and the number of cells with `NA` values

## See also

[`Summary-methods`](https://rspatial.github.io/terra/reference/summarize-generics.md),
[`global`](https://rspatial.github.io/terra/reference/global.md),
[`quantile`](https://rspatial.github.io/terra/reference/quantile.md)

## Examples

``` r
set.seed(0)
r <- rast(nrows=10, ncols=10, nlyrs=3)
values(r) <- runif(nlyr(r)*ncell(r))
summary(r)
#>      lyr.1             lyr.2             lyr.3        
#>  Min.   :0.01339   Min.   :0.01308   Min.   :0.02779  
#>  1st Qu.:0.32308   1st Qu.:0.28440   1st Qu.:0.18947  
#>  Median :0.48781   Median :0.51860   Median :0.37810  
#>  Mean   :0.52076   Mean   :0.51569   Mean   :0.43403  
#>  3rd Qu.:0.77171   3rd Qu.:0.72570   3rd Qu.:0.63596  
#>  Max.   :0.99191   Max.   :0.99268   Max.   :0.98156  
```
