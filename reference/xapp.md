# Apply a function to the cells of two SpatRasters

Apply a function to the values of each cell of two (multilayer)
SpatRasters.

## Usage

``` r
# S4 method for class 'SpatRaster,SpatRaster'
xapp(x, y, fun, ..., filename="", overwrite=FALSE, wopt=list())
```

## Arguments

- x:

  SpatRaster

- y:

  SpatRaster with the same geometry as `x`

- fun:

  a function that operates on two vectors

- ...:

  additional arguments for `fun`. These are typically numerical
  constants. They should \*never\* be another SpatRaster

- filename:

  character. Output filename

- overwrite:

  logical. If `TRUE`, `filename` is overwritten

- wopt:

  list with named options for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## Value

SpatRaster

## See also

[`app`](https://rspatial.github.io/terra/reference/app.md),
[`lapp`](https://rspatial.github.io/terra/reference/lapp.md),
[`tapp`](https://rspatial.github.io/terra/reference/tapp.md),
[`Math-methods`](https://rspatial.github.io/terra/reference/math-generics.md),
[`roll`](https://rspatial.github.io/terra/reference/roll.md)

## Examples

``` r
r <- rast(ncols=10, nrows=10, nlyr=5)
set.seed(1)
r <- init(r, runif)
s <- init(r, runif)
x <- xapp(r, s, fun=cor)
```
