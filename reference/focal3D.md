# Three-dimensional focal values

Calculate focal ("moving window") values for the three-dimensional
neighborhood (window) of focal cells. See
[`focal`](https://rspatial.github.io/terra/reference/focal.md) for
two-dimensional focal computation.

## Usage

``` r
# S4 method for class 'SpatRaster'
focal3D(x, w=3, fun=mean, ..., na.policy="all", fillvalue=NA, pad=FALSE, 
  padvalue=fillvalue, expand=FALSE, silent=TRUE, 
  filename="", overwrite=FALSE, wopt=list())
```

## Arguments

- x:

  SpatRaster

- w:

  window. A rectangular prism (cuboid) defined by three numbers or by a
  three-dimensional array. The values are used as weights, and are
  usually zero, one, NA, or fractions. The window used must have odd
  dimensions. If you desire to use even sides, you can use an array, and
  pad the values with rows and/or columns that contain only `NA`s.

- fun:

  function that takes multiple numbers, and returns one or multiple
  numbers for each focal area. For example mean, modal, min or max

- ...:

  additional arguments passed to `fun` such as `na.rm`

- na.policy:

  character. Can be used to determine the cells of `x`, in the central
  layer, for which focal values should be computed. Must be one of "all"
  (compute for all cells), "only" (only for cells that are `NA`) or
  "omit" (skip cells that are `NA`). Note that the value of this
  argument does not affect which cells around each focal cell are
  included in the computations (use `na.rm=TRUE` to ignore cells that
  are `NA` in the computation of the focal value)

- fillvalue:

  numeric. The value of the cells in the virtual rows and columns
  outside of the raster

- pad:

  logical. Add virtual layers before the first and after the last layer

- padvalue:

  numeric. The value of the cells in the virtual layers

- expand:

  logical. Add virtual layers before the first or after the last layer
  that are the same as the first or last layers. If `TRUE`, arguments
  `pad` and `padvalue` are ignored

- silent:

  logical. If `TRUE` error messages are printed that may occur when
  trying `fun` to determine the length of the returned value. This can
  be useful in debugging a function passed to `fun` that does not work

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

[`focal`](https://rspatial.github.io/terra/reference/focal.md)

## Examples

``` r
r <- rast(system.file("ex/logo.tif", package="terra"))   
x <- focal3D(r, c(5,5,3), na.rm=TRUE)

a <- array(c(0,1,0,1,1,1,0,1,0, rep(1,9), 0,1,0,1,1,1,0,1,0), c(3,3,3))
a[a==0] <- NA
z <- focal3D(r, a, na.rm=TRUE)
```
