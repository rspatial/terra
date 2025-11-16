# Initialize a SpatRaster with values

Create a SpatRaster with values reflecting a cell property: "x", "y",
"xy", "col", "row", "cell" or "chess". Alternatively, a function can be
used. In that case, cell values are initialized without reference to
pre-existing values. E.g., initialize with a random number
(`fun=`[`runif`](https://rdrr.io/r/stats/Uniform.html)). While there are
more direct ways of achieving this for small objects (see examples) for
which a vector with all values can be created in memory, the `init`
function will also work for SpatRasters with many cells.

## Usage

``` r
# S4 method for class 'SpatRaster'
init(x, fun, ..., filename="", overwrite=FALSE, wopt=list())
```

## Arguments

- x:

  SpatRaster

- fun:

  function to be applied. This must be a either single number, multiple
  numbers, a function, or one of a set of known character values. A
  function must take the number of cells as a single argument to return
  a vector of values with a length equal to the number of cells, such as
  `fun=runif`. Allowed character values are "x", "y", "row", "col",
  "cell", and "chess" to get the x or y coordinate or both, row, col or
  cell number or a chessboard pattern (alternating 0 and 1 values)

- ...:

  additional arguments passed to `fun`

- filename:

  character. Output filename

- overwrite:

  logical. If `TRUE`, `filename` is overwritten

- wopt:

  additional arguments for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## Value

SpatRaster

## Examples

``` r
r <- rast(ncols=10, nrows=5, xmin=0, xmax=10, ymin=0, ymax=5)
x <- init(r, fun="cell")
y <- init(r, fun=runif)

# initialize with a single value 
z <- init(r, fun=8)
```
