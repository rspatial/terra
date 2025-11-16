# Focal regression

Calculate values for a moving-window by comparing the value in one
layers with the values in one to many other layers. A typical case is
the computation of the coefficients for a focal linear regression model.

## Usage

``` r
# S4 method for class 'SpatRaster'
focalReg(x, w=3, fun="ols", ..., fillvalue=NA, filename="", overwrite=FALSE, wopt=list())
```

## Arguments

- x:

  SpatRaster with at least two layers. The first is the "Y" (dependent)
  variable and the remainder are the "X" (independent) variables

- w:

  numeric or matrix to define the focal window. The window an be defined
  as one (for a square) or two numbers (row, col); or with an odd-sized
  weights matrix. See the Details section in
  [`focal`](https://rspatial.github.io/terra/reference/focal.md). Note
  that if a matrix with numbers other than zero or one are used, the
  values are used as weights. For this to work, `fun` must have an
  argument `weights`

- fun:

  a function with at least two arguments (one for each layer). There is
  a built-in function "ols" for both the weighted and unweighted
  Ordinary Least Square regression. This function has an additional
  argument `na.rm=FALSE` and `intercept=TRUE`

- ...:

  additional arguments for `fun`

- fillvalue:

  numeric. The value of the cells in the virtual rows and columns
  outside of the raster

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

[`focal`](https://rspatial.github.io/terra/reference/focal.md)`, `[`focal3D`](https://rspatial.github.io/terra/reference/focal3D.md),
[focalValues](https://rspatial.github.io/terra/reference/focalValues.md)

## Examples

``` r
r <- rast(ncols=10, nrows=10, ext(0, 10, 0, 10))
values(r) <- 1:ncell(r)
x <- c(r, init(r, runif) * r)
f <- focalReg(x, 3)
```
