# Origin

Get or set the coordinates of the point of origin of a SpatRaster. This
is the point closest to (0, 0) that you could get if you moved towards
that point in steps of the x and y resolution.

## Usage

``` r
# S4 method for class 'SpatRaster'
origin(x)

# S4 method for class 'SpatRaster'
origin(x) <- value
```

## Arguments

- x:

  SpatRaster

- value:

  numeric vector of length 1 or 2

## Value

A vector of two numbers (x and y coordinates)

## Examples

``` r
r <- rast(xmin=-0.5, xmax = 9.5, ncols=10)
origin(r)
#> [1] 0.5 0.0
origin(r) <- c(0,0)
r
#> class       : SpatRaster 
#> size        : 180, 10, 1  (nrow, ncol, nlyr)
#> resolution  : 1, 1  (x, y)
#> extent      : 0, 10, -90, 90  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (CRS84) (OGC:CRS84) 
```
