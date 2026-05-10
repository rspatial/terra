# Create a text representation of (the skeleton of) an object

Create a text representation of (the skeleton of) an object

## Usage

``` r
# S4 method for class 'SpatExtent'
as.character(x)

# S4 method for class 'SpatRaster'
as.character(x)
```

## Arguments

- x:

  SpatRaster or SpatExtent

## Value

character

## Examples

``` r
r <- rast()
r
#> class       : SpatRaster
#> size        : 180, 360, 1  (nrow, ncol, nlyr)
#> resolution  : 1, 1  (x, y)
#> extent      : -180, 180, -90, 90  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (CRS84) (OGC:CRS84)
as.character(r)
#> [1] "rast(ncols=360, nrows=180, nlyrs=1, xmin=-180, xmax=180, ymin=-90, ymax=90, names=c('lyr.1'), crs='OGC:CRS84')"

e <- ext(r)
e
#> SpatExtent : -180, 180, -90, 90 (xmin, xmax, ymin, ymax)
as.character(e)
#> [1] "ext(-180, 180, -90, 90)"
```
