# Trim a SpatRaster

Trim (shrink) a SpatRaster by removing outer rows and columns that are
`NA` or another value.

## Usage

``` r
# S4 method for class 'SpatRaster'
trim(x, padding=0, value=NA, filename="", ...)
```

## Arguments

- x:

  SpatRaster

- padding:

  integer. Number of outer rows/columns to keep

- value:

  numeric. The value of outer rows or columns that are to be removed

- filename:

  character. Output filename

- ...:

  additional arguments for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## Value

SpatRaster

## See also

[`extend`](https://rspatial.github.io/terra/reference/extend.md)

## Examples

``` r
r <- rast(ncols=10, nrows=10, xmin=0,xmax=10,ymin=0,ymax=10)
v <- rep(NA, ncell(r))
v[c(12,34,69)] <- 1:3
values(r) <- v
s <- trim(r) 
```
