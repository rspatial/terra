# is not NA

Shortcut method to avoid the two-step `!is.na(x)`

## Usage

``` r
# S4 method for class 'SpatRaster'
not.na(x, falseNA=FALSE, filename="", ...)
```

## Arguments

- x:

  SpatRaster

- falseNA:

  logical. If `TRUE`, the output cell values are either `TRUE`, for
  cells that are not `NA` in `x`, or `NA` for the cells that are `NA` in
  `x`. Otherwise, the output values are either `TRUE` or `FALSE`

- filename:

  character. Output filename

- ...:

  additional arguments for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## See also

[`Compare-methods`](https://rspatial.github.io/terra/reference/compare-generics.md)

## Value

SpatRaster

## Examples

``` r
r <- rast(ncols=5, nrows=5, vals=1, ext=c(0,1,0,1))
r[10:20] <- NA
x <- not.na(r)
y <- not.na(r, falseNA=TRUE)
unique(values(c(x, y)))
#>      lyr.1 lyr.1
#> [1,]  TRUE  TRUE
#> [2,] FALSE    NA
```
