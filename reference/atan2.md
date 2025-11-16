# Two argument arc-tangent

For SpatRasters x and y, atan2(y, x) returns the angle in radians for
the tangent y/x, handling the case when x is zero. See
[`Trig`](https://rdrr.io/r/base/Trig.html)

See
[`Math-methods`](https://rspatial.github.io/terra/reference/math-generics.md)
for other trigonometric and mathematical functions that can be used with
SpatRasters.

## Usage

``` r
# S4 method for class 'SpatRaster,SpatRaster'
atan2(y, x)

# S4 method for class 'SpatRaster,SpatRaster'
atan_2(y, x, filename, ...)
```

## Arguments

- y:

  SpatRaster

- x:

  SpatRaster

- filename:

  character. Output filename

- ...:

  additional arguments for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## See also

[`Math-methods`](https://rspatial.github.io/terra/reference/math-generics.md)

## Examples

``` r
r1 <- rast(nrows=10, ncols=10)
r2 <- rast(nrows=10, ncols=10)
values(r1) <- (runif(ncell(r1))-0.5) * 10
values(r2) <- (runif(ncell(r1))-0.5) * 10
atan2(r1, r2)
#> class       : SpatRaster 
#> size        : 10, 10, 1  (nrow, ncol, nlyr)
#> resolution  : 36, 18  (x, y)
#> extent      : -180, 180, -90, 90  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (CRS84) (OGC:CRS84) 
#> source(s)   : memory
#> name        :     lyr.1 
#> min value   : -3.125421 
#> max value   :  3.077551 
```
