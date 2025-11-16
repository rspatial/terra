# Add (in place) a SpatRaster to another SpatRaster or to a SpatRasterDataset or SpatRasterCollection

Add (in place) a SpatRaster to another SpatRaster. Comparable with
[`c`](https://rspatial.github.io/terra/reference/c.md), but without
copying the object.

## Usage

``` r
# S4 method for class 'SpatRaster,SpatRaster'
add(x) <- value

# S4 method for class 'SpatRasterDataset,SpatRaster'
add(x) <- value

# S4 method for class 'SpatRasterCollection,SpatRaster'
add(x) <- value
```

## Arguments

- x:

  SpatRaster, SpatRasterDataset or SpatRasterCollection

- value:

  SpatRaster

## See also

[`c`](https://rspatial.github.io/terra/reference/c.md)

## Value

SpatRaster

## Examples

``` r
r <- rast(nrows=5, ncols=9, vals=1:45)
x <- c(r, r*2)
add(x) <- r*3
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
```
