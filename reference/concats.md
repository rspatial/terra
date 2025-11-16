# Concatenate categorical rasters

Combine two categorical rasters by concatenating their levels.

## Usage

``` r
# S4 method for class 'SpatRaster'
concats(x, y, filename="", ...)
```

## Arguments

- x:

  SpatRaster (with a single, categorical, layer)

- y:

  SpatRaster (with a single, categorical, layer)

- filename:

  character. Output filename

- ...:

  additional arguments for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## Value

SpatRaster

## See also

[`cats`](https://rspatial.github.io/terra/reference/factors.md)

## Examples

``` r
set.seed(0)
r <- rast(nrows=10, ncols=10)
values(r) <- sample(3, ncell(r), replace=TRUE)
levels(r) <- data.frame(id=1:3, cover=c("forest", "water", "urban"))

rr <- rast(r)
values(rr) <- sample(1:3, ncell(rr), replace=TRUE)
levels(rr) <- data.frame(id=c(1:3), color=c("red", "green", "blue"))

x <- concats(r, rr)
x
#> class       : SpatRaster 
#> size        : 10, 10, 1  (nrow, ncol, nlyr)
#> resolution  : 36, 18  (x, y)
#> extent      : -180, 180, -90, 90  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (CRS84) (OGC:CRS84) 
#> source(s)   : memory
#> categories  : cover_color 
#> name        : cover_color 
#> min value   :  forest_red 
#> max value   :  urban_blue 
levels(x)[[1]]
#>   ID  cover_color
#> 1  0   forest_red
#> 2  1 forest_green
#> 3  2  forest_blue
#> 4  3    water_red
#> 5  4  water_green
#> 6  5   water_blue
#> 7  6    urban_red
#> 8  7  urban_green
#> 9  8   urban_blue
```
