# Set a window

Assign a window (area of interest) to a SpatRaster with a `SpatExtent`,
or set it to `NULL` to remove the window. This is similar to
[`crop`](https://rspatial.github.io/terra/reference/crop.md) without
actually creating a new dataset.

The window is intersect with the extent of the SpatRaster. It is
envisioned that in a future version, the window may also go outside
these boundaries.

## Usage

``` r
# S4 method for class 'SpatRaster'
window(x) <- value

# S4 method for class 'SpatRaster'
window(x)
```

## Arguments

- x:

  SpatRaster

- value:

  SpatExtent

## See also

[`crop`](https://rspatial.github.io/terra/reference/crop.md),
[`extend`](https://rspatial.github.io/terra/reference/extend.md)

## Value

none for `window<-` and logical for `window`

## Examples

``` r
f <- system.file("ex/elev.tif", package="terra") 
r <- rast(f)
global(r, "mean", na.rm=TRUE)
#>               mean
#> elevation 348.3366
e <- ext(c(5.9, 6,49.95, 50))

window(r) <- e
global(r, "mean", na.rm=TRUE)
#>               mean
#> elevation 396.4167
r
#> class       : SpatRaster 
#> size        : 6, 12, 1  (nrow, ncol, nlyr)
#> resolution  : 0.008333333, 0.008333333  (x, y)
#> window      : 5.9, 6, 49.95, 50  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (EPSG:4326) 
#> source      : elev.tif 
#> name        : elevation 
#> min value   :      >141 
#> max value   :      547< 

x <- rast(f)
xe <- crop(x, e)
global(xe, "mean", na.rm=TRUE)
#>               mean
#> elevation 396.4167

b <- c(xe, r)
window(b)
#> [1] FALSE  TRUE
b
#> class       : SpatRaster 
#> size        : 6, 12, 2  (nrow, ncol, nlyr)
#> resolution  : 0.008333333, 0.008333333  (x, y)
#> extent (win): 5.9, 6, 49.95, 50  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (EPSG:4326) 
#> sources     : memory  
#>               elev.tif  
#> varnames    : elev 
#>               elev 
#> names       : elevation, elevation 
#> min values  :       323,      >141 
#> max values  :       492,      547< 

window(r) <- NULL
r
#> class       : SpatRaster 
#> size        : 90, 95, 1  (nrow, ncol, nlyr)
#> resolution  : 0.008333333, 0.008333333  (x, y)
#> extent      : 5.741667, 6.533333, 49.44167, 50.19167  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (EPSG:4326) 
#> source      : elev.tif 
#> name        : elevation 
#> min value   :       141 
#> max value   :       547 
```
