# global statistics

Compute global statistics, that is summarized values of an entire
SpatRaster.

If `x` is very large `global` can fail, except when `fun` is one of
these built-in functions "mean", "min", "max", "sum", "prod", "range"
(min and max), "rms" (root mean square), "sd" (sample standard
deviation), "std" (population standard deviation), "isNA" (number of
cells that are NA), "notNA" (number of cells that are not NA), "anyNA",
"anynotNA". Note that "anyNA" and "anynotNA" cannot be combined with
other functions.

The reason that this can fail with large raster and a custom function is
that all values need to be loaded into memory. To circumvent this
problem you can run `global` with a sample of the cells.

You can compute a weighted mean or sum by providing a SpatRaster with
weights.

## Usage

``` r
# S4 method for class 'SpatRaster'
global(x, fun="mean", weights=NULL, maxcell=Inf, ...)
```

## Arguments

- x:

  SpatRaster

- fun:

  function to be applied to summarize the values by zone. Either as one
  or more of these built-in character values: "max", "min", "mean",
  "sum", "range", "rms" (root mean square), "sd", "std" (population sd,
  using `n` rather than `n-1`), "isNA", "notNA", "anyNA", "anynotNA"; or
  a proper R function (but these may fail for very large SpatRasters
  unless you specify `maxcell`)

- ...:

  additional arguments passed on to `fun`

- weights:

  NULL or SpatRaster

- maxcell:

  positive integer used to take a regular sample of `x`. Ignored by the
  built-in functions.

## Value

A `data.frame` with a row for each layer

## See also

[`zonal`](https://rspatial.github.io/terra/reference/zonal.md) for
"zonal" statistics, and
[`app`](https://rspatial.github.io/terra/reference/app.md) or
[`Summary-methods`](https://rspatial.github.io/terra/reference/summarize-generics.md)
for "local" statistics, and
[`extract`](https://rspatial.github.io/terra/reference/extract.md) for
summarizing values for polygons. Also see
[`focal`](https://rspatial.github.io/terra/reference/focal.md) for
"focal" or "moving window" operations.

## Examples

``` r
r <- rast(ncols=10, nrows=10)
values(r) <- 1:ncell(r)
global(r, "sum")
#>        sum
#> lyr.1 5050
global(r, "mean", na.rm=TRUE)
#>       mean
#> lyr.1 50.5
x <- c(r, r/10)
global(x, c("sum", "mean", "sd"), na.rm=TRUE)
#>          sum  mean        sd
#> lyr.1   5050 50.50 29.011492
#> lyr.1.1  505  5.05  2.901149

global(x, function(i) min(i) / max(i))
#>         global
#> lyr.1     0.01
#> lyr.1.1   0.01
```
