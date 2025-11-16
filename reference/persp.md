# Perspective plot

Perspective plot of a SpatRaster. This is an implementation of a generic
function in the graphics package.

## Usage

``` r
# S4 method for class 'SpatRaster'
persp(x, maxcells=100000, ...)
```

## Arguments

- x:

  SpatRaster. Only the first layer is used

- maxcells:

  integer \> 0. Maximum number of cells to use for the plot. If
  `maxpixels < ncell(x)`, `spatSample(method="regular")` is used before
  plotting

- ...:

  Any argument that can be passed to
  [`persp`](https://rdrr.io/r/graphics/persp.html) (graphics package)

## See also

[`persp`](https://rdrr.io/r/graphics/persp.html), `contour`, `plot`

## Examples

``` r
r <- rast(system.file("ex/elev.tif", package="terra"))
persp(r)
```
