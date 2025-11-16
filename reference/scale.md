# Scale values

Center and/or scale raster data. For details see
[`scale`](https://rdrr.io/r/base/scale.html)

## Usage

``` r
# S4 method for class 'SpatRaster'
scale(x, center=TRUE, scale=TRUE)
```

## Arguments

- x:

  SpatRaster

- center:

  logical or numeric. If `TRUE`, centering is done by subtracting the
  layer means (omitting `NA`s), and if `FALSE`, no centering is done. If
  `center` is a numeric vector (recycled to `nlyr(x)`), then each layer
  of `x` has the corresponding value from center subtracted from it.

- scale:

  logical or numeric. If `TRUE`, scaling is done by dividing the
  (centered) layers of `x` by their standard deviations if `center` is
  `TRUE`, and the root mean square otherwise. If scale is `FALSE`, no
  scaling is done. If `scale` is a numeric vector (recycled to
  `nlyr(x)`), each layer of `x` is divided by the corresponding value.
  Scaling is done after centering.

## Value

SpatRaster

## See also

[`scale_linear`](https://rspatial.github.io/terra/reference/scale_linear.md)

## Examples

``` r
r <- rast(system.file("ex/logo.tif", package="terra"))   
s <- scale(r)

## the equivalent, computed in steps
m <- global(r, "mean")
rr <- r - m[,1]
rms <- global(rr, "rms")
ss <- rr / rms[,1]
```
