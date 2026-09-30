# Distance to reference values

For each cell, compute the distance in value-space to reference values.
References are a matrix or data.frame whose columns match the layer
names of `x`, or a SpatVector of locations from which values are
extracted.

## Usage

``` r
# S4 method for class 'SpatRaster,matrix'
distValues(x, y, fun="squared", weights=NULL, ...,
    filename="", overwrite=FALSE, wopt=list())

# S4 method for class 'SpatRaster,data.frame'
distValues(x, y, fun="squared", weights=NULL, ...,
    filename="", overwrite=FALSE, wopt=list())

# S4 method for class 'SpatRaster,SpatVector'
distValues(x, y, fun="squared", center=FALSE, scale=FALSE,
    weights=NULL, ..., filename="", overwrite=FALSE, wopt=list())
```

## Arguments

- x:

  SpatRaster

- y:

  matrix, data.frame, or SpatVector

- fun:

  character. `"abs"` (mean absolute difference) or `"squared"` (mean
  squared difference). Or a function

- weights:

  numeric. Optional weights for the layers of `x`

- center:

  logical or numeric. Passed to
  [`scale`](https://rspatial.github.io/terra/reference/scale.md)
  (SpatVector method only)

- scale:

  logical or numeric. Passed to
  [`scale`](https://rspatial.github.io/terra/reference/scale.md)
  (SpatVector method only)

- ...:

  additional arguments passed to `fun` (for the built-in functions,
  `na.rm=TRUE`)

- filename:

  character. Output filename

- overwrite:

  logical. If `TRUE`, `filename` is overwritten

- wopt:

  additional arguments for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## Value

SpatRaster

## See also

[`bestMatch`](https://rspatial.github.io/terra/reference/bestMatch.md),
[`scale`](https://rspatial.github.io/terra/reference/scale.md)

## Examples

``` r
r <- rast(system.file("ex/logo.tif", package="terra"))
pts <- vect(cbind(c(25.25, 34.324), c(54.577, 46.489)))
x <- scale(r)
d <- distValues(x, pts)
```
