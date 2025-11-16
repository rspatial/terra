# Where are the cells with the min or max values?

This method returns the cell numbers for the cells with the min or max
values of each layer in a SpatRaster.

## Usage

``` r
# S4 method for class 'SpatRaster'
where.min(x, values=TRUE, list=FALSE)

# S4 method for class 'SpatRaster'
where.max(x, values=TRUE, list=FALSE)
```

## Arguments

- x:

  SpatRaster

- values:

  logical. If `TRUE` the min or max values are also returned

- list:

  logical. If `TRUE` a list is returned instead of a matrix

## Value

matrix or list

## See also

[`which`](https://rdrr.io/r/base/which.html) and
[`Summary-methods`](https://rspatial.github.io/terra/reference/summarize-generics.md)
for `which.min` and `which.max`

## Examples

``` r
f <- system.file("ex/elev.tif", package="terra")
r <- rast(f)
where.min(r)
#>      layer cell value
#> [1,]     1 7770   141
#> [2,]     1 8055   141
```
