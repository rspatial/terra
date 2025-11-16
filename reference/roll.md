# Rolling (moving) functions

Compute "rolling" or "moving" values, such as the "rolling average" for
each cell in a SpatRaster.

See [`focal`](https://rspatial.github.io/terra/reference/focal.md) for
spatially moving averages and similar computations. And see
[`cumsum`](https://rspatial.github.io/terra/reference/math-generics.md)
and other cum\* functions to compute cumulate values.

## Usage

``` r
# S4 method for class 'SpatRaster'
roll(x, n, fun=mean, type="around", circular=FALSE,
      na.rm=FALSE, filename="", ..., wopt=list()) 

# S4 method for class 'numeric'
roll(x, n, fun=mean, type="around", circular=FALSE, na.rm=FALSE, ...)
```

## Arguments

- x:

  SpatRaster or numeric

- n:

  integer \> 1. The size of the "window", that is, the number of
  sequential cells to use in `fun`

- fun:

  a function like mean, min, max, sum

- type:

  character. One of "around", "to", or "from". The choice indicates
  which values should be used in the computation. The focal cell is
  always used. If type is "around", `(n-1)/2` before and after the focal
  cell are also included. If type = "from", `n-1` cells are after the
  focal cell are included. If type = "to", `n-1` cells before the focal
  cell are included. For example, when using n=3 for element 5 of a
  vector; "around" used elements 4,5,6; "to" used elements 3,4,5, and
  "from" uses elements 5,6,7

- circular:

  logical. If `TRUE`, the data are considered to have a circular nature
  (e.g. days or months of the year), such that there are no missing
  values before first or after the last value.

- na.rm:

  logical. If `TRUE`, `NA` values should be ignored (by `fun`)

- filename:

  character. Output filename

- ...:

  additional arguments for `fun`

- wopt:

  list with named options for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## Value

Same as `x`

## See also

[`cumsum`](https://rspatial.github.io/terra/reference/math-generics.md),
[`focal`](https://rspatial.github.io/terra/reference/focal.md)

## Examples

``` r
## numeric
roll(1:12, 3, mean)
#>  [1] NA  2  3  4  5  6  7  8  9 10 11 NA
roll(1:12, 3, mean, "to")
#>  [1] NA NA  2  3  4  5  6  7  8  9 10 11
roll(1:12, 3, mean, circular=TRUE)
#>  [1]  5  2  3  4  5  6  7  8  9 10 11  8

## SpatRaster
r <- rast(ncol=2, nrow=2, nlyr=10, vals=1)
r[1,2] <- 2
r[2,2] <- 4

values(roll(r, n=3, "sum", "from", na.rm=FALSE))
#>      lyr.1 lyr.2 lyr.3 lyr.4 lyr.5 lyr.6 lyr.7 lyr.8 lyr.9 lyr.10
#> [1,]     3     3     3     3     3     3     3     3   NaN    NaN
#> [2,]     6     6     6     6     6     6     6     6   NaN    NaN
#> [3,]     3     3     3     3     3     3     3     3   NaN    NaN
#> [4,]    12    12    12    12    12    12    12    12   NaN    NaN
values(roll(r, n=3, "sum", "from", na.rm=TRUE))
#>      lyr.1 lyr.2 lyr.3 lyr.4 lyr.5 lyr.6 lyr.7 lyr.8 lyr.9 lyr.10
#> [1,]     3     3     3     3     3     3     3     3     2      1
#> [2,]     6     6     6     6     6     6     6     6     4      2
#> [3,]     3     3     3     3     3     3     3     3     2      1
#> [4,]    12    12    12    12    12    12    12    12     8      4
values(roll(r, n=3, "sum", "from", circular=TRUE))
#>      lyr.1 lyr.2 lyr.3 lyr.4 lyr.5 lyr.6 lyr.7 lyr.8 lyr.9 lyr.10
#> [1,]     3     3     3     3     3     3     3     3     3      3
#> [2,]     6     6     6     6     6     6     6     6     6      6
#> [3,]     3     3     3     3     3     3     3     3     3      3
#> [4,]    12    12    12    12    12    12    12    12    12     12

values(roll(r, n=3, "sum", "to", na.rm=TRUE))
#>      lyr.1 lyr.2 lyr.3 lyr.4 lyr.5 lyr.6 lyr.7 lyr.8 lyr.9 lyr.10
#> [1,]     1     2     3     3     3     3     3     3     3      3
#> [2,]     2     4     6     6     6     6     6     6     6      6
#> [3,]     1     2     3     3     3     3     3     3     3      3
#> [4,]     4     8    12    12    12    12    12    12    12     12

values(roll(r, n=3, "sum", "around", circular=TRUE))
#>      lyr.1 lyr.2 lyr.3 lyr.4 lyr.5 lyr.6 lyr.7 lyr.8 lyr.9 lyr.10
#> [1,]     3     3     3     3     3     3     3     3     3      3
#> [2,]     6     6     6     6     6     6     6     6     6      6
#> [3,]     3     3     3     3     3     3     3     3     3      3
#> [4,]    12    12    12    12    12    12    12    12    12     12
```
