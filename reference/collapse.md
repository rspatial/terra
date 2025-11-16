# tighten SpatRaster or SpatRasterDataset objects

Combines data sources within a SpatRaster (that are in memory, or from
the same file) to allow for faster processing.

Or combine sub-datasets into a SpatRaster.

## Usage

``` r
# S4 method for class 'SpatRaster'
tighten(x)

# S4 method for class 'SpatRasterDataset'
tighten(x)
```

## Arguments

- x:

  SpatRaster or SpatRasterDataset

## Value

SpatRaster

## Examples

``` r
r <- rast(nrow=5, ncol=9, vals=1:45)
x <- c(r, r*2, r*3)
x
#> class       : SpatRaster 
#> size        : 5, 9, 3  (nrow, ncol, nlyr)
#> resolution  : 40, 36  (x, y)
#> extent      : -180, 180, -90, 90  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (CRS84) (OGC:CRS84) 
#> source(s)   : memory
#> names       : lyr.1, lyr.1, lyr.1 
#> min values  :     1,     2,     3 
#> max values  :    45,    90,   135 
tighten(x)
#> class       : SpatRaster 
#> size        : 5, 9, 3  (nrow, ncol, nlyr)
#> resolution  : 40, 36  (x, y)
#> extent      : -180, 180, -90, 90  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (CRS84) (OGC:CRS84) 
#> source(s)   : memory
#> names       : lyr.1, lyr.1, lyr.1 
#> min values  :     1,     2,     3 
#> max values  :    45,    90,   135 
```
