# Make a SpatRaster method memory-safe

This method allows for running a function that takes a SpatRaster as
first argument in chunks (tiles). This can be useful if the functions is
not memory-safe, typically because it reads all the raster cell values
into memory.

This method is not designed to be especially efficient, and there might
be more efficient ways to accomplish what the the goal of the function
that is not memory-safe.

Also, some functions must have access to all cells at once to be valid.
In those cases, chunk would return incorrect results.

## Usage

``` r
# S4 method for class 'SpatRaster'
chunk(x, fun, ..., n=NULL, buffer=0, filename="", wopt=list())
```

## Arguments

- x:

  SpatRaster

- fun:

  function that takes a SpatRaster as first argument

- ...:

  additional arguments for fun

- n:

  NULL or positive integer to specifying the number of rows and columns
  for each chunk (or 2 numbers for a different number of rows and
  columns, as in
  [`getTileExtents`](https://rspatial.github.io/terra/reference/makeTiles.md)

- buffer:

  integer. The number of additional rows and columns added to each tile.
  Can be a single number, or two numbers to specify a separate number of
  rows and columns. This allows for creating overlapping tiles that can
  be used for computing spatial context dependent values with e.g.
  [`focal`](https://rspatial.github.io/terra/reference/focal.md). The
  expansion is only inside `x`, no rows or columns outside of `x` are
  added

- filename:

  character. Output filename

- wopt:

  list with additional arguments for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## Value

SpatRaster

## Examples

``` r
s <- rast(system.file("ex/logo.tif", package="terra"))
f <- function(x, a = 0) {
  print("chunk")
  sum(x) + a
}

x <- chunk(s, f, a=100)
#> [1] "chunk"
#> [1] "chunk"
#> [1] "chunk"
#> [1] "chunk"
#> [1] "chunk"
#> [1] "chunk"
```
