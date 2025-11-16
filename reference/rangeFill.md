# Fill layers with a range

Fill layers with cell-varying ranges defined by a start and end
SpatRaster. The range must start at 1 and end at a user-defined maximum.
Output values are either zero (not in the range) or one (in the range).

For example, for a cell with `start=3`, `end=5` and with `limit=8`, the
output for that cell would be `0,0,1,1,1,0,0,0`

## Usage

``` r
# S4 method for class 'SpatRaster'
rangeFill(x, limit, circular=FALSE, filename="", ...)
```

## Arguments

- x:

  SpatRaster with at two layers. The cell values of the first layer
  indicate the start of the range (1 based); the cell values are
  indicate the end of the range

- limit:

  numeric \> 1. The range size

- circular:

  logical. If `TRUE` the values are considered circular, such as the
  days of the year. In that case, if first \> last the layers used are
  c(first:limit, 1:last). Otherwise, if `circular=FALSE`, such a range
  would be considered invalid and `NA` would be used

- filename:

  character. Output filename

- ...:

  additional arguments for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## Value

SpatRaster

## See also

[`rapp`](https://rspatial.github.io/terra/reference/rapp.md)

## Examples

``` r
x <- y <- rast(ncol=2, nrow=2)
values(x) <- c(NA, 1:3)
values(y) <- c(NA, 4:6)

r <- rangeFill(c(x, y), 8)
```
