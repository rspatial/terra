# depth of SpatRaster layers

Get or set the depth of the layers of a SpatRaster. Experimental.

## Usage

``` r
# S4 method for class 'SpatRaster'
depth(x)

# S4 method for class 'SpatRaster'
depth(x) <- value

# S4 method for class 'SpatRaster'
depthName(x)

# S4 method for class 'SpatRaster'
depthName(x) <- value

# S4 method for class 'SpatRaster'
depthUnit(x)

# S4 method for class 'SpatRaster'
depthUnit(x) <- value
```

## Arguments

- x:

  SpatRaster

- value:

  numeric vector (`depth`), or character (`depthName` and `depthUnit`)

## Value

numeric or SpatRaster

## See also

[`time`](https://rspatial.github.io/terra/reference/time.md)

## Examples

``` r
s <- rast(nlyr=3)

depth(s) <- c(0, pi/2, pi)
depth(s)
#> [1] 0.000000 1.570796 3.141593

depthName(s) <- "angle"
depthUnit(s) <- "radians"
s
#> class       : SpatRaster 
#> size        : 180, 360, 3  (nrow, ncol, nlyr)
#> resolution  : 1, 1  (x, y)
#> extent      : -180, 180, -90, 90  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (CRS84) (OGC:CRS84) 
#> depth       : 0 to 3.14159265358979 (angle [radians]: 3 steps) 
```
