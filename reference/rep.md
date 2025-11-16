# Replicate layers

Replicate layers in a SpatRaster

## Usage

``` r
# S4 method for class 'SpatRaster'
rep(x, ...)
```

## Arguments

- x:

  SpatRaster

- ...:

  arguments as in [`rep`](https://rdrr.io/r/base/rep.html)

## Value

SpatRaster

## Examples

``` r
s <- rast(system.file("ex/logo.tif", package="terra"))   
x <- rep(s, 2)
nlyr(x)
#> [1] 6
names(x)
#> [1] "red"   "green" "blue"  "red"   "green" "blue" 
x
#> class       : SpatRaster 
#> size        : 77, 101, 6  (nrow, ncol, nlyr)
#> resolution  : 1, 1  (x, y)
#> extent      : 0, 101, 0, 77  (xmin, xmax, ymin, ymax)
#> coord. ref. : Cartesian (Meter) 
#> source      : logo.tif 
#> names       : red, green, blue, red, green, blue 
#> min values  :   0,     0,    0,   0,     0,    0 
#> max values  : 255,   255,  255, 255,   255,  255 
```
