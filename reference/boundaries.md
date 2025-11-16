# Detect boundaries (edges)

Detect boundaries (edges). Boundaries are cells that have more than one
class in the 4 or 8 cells surrounding it, or, if `classes=FALSE`, cells
with values and cells with `NA`.

## Usage

``` r
# S4 method for class 'SpatRaster'
boundaries(x, classes=FALSE, inner=TRUE, 
         directions=8, falseval=0, filename="", ...)
```

## Arguments

- x:

  SpatRaster

- inner:

  logical. If `TRUE`, "inner" boundaries are returned, else "outer"
  boundaries are returned

- classes:

  character. Logical. If `TRUE` all different values are (after
  rounding) distinguished, as well as `NA`. If `FALSE` (the default)
  only edges between `NA` and non-`NA` cells are considered

- directions:

  integer. Which cells are considered adjacent? Should be 8 (Queen's
  case) or 4 (Rook's case)

- falseval:

  numeric. The value to use for cells that are not a boundary and not
  `NA`

- filename:

  character. Output filename

- ...:

  options for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## Value

SpatRaster. Cell values are either 1 (a border) or 0 (not a border), or
`NA`

## See also

[`focal`](https://rspatial.github.io/terra/reference/focal.md),
[`patches`](https://rspatial.github.io/terra/reference/patches.md)

## Examples

``` r
r <- rast(nrows=18, ncols=36, xmin=0)
r[150:250] <- 1
r[251:450] <- 2
bi <- boundaries(r)
bo <- boundaries(r, inner=FALSE)
bc <- boundaries(r, classes=TRUE)
#plot(bc)
```
