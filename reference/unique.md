# Unique values

This method returns the unique values in a SpatRaster, or removes
duplicates records (geometry and attributes) in a SpatVector.

## Usage

``` r
# S4 method for class 'SpatRaster'
unique(x, incomparables=FALSE, digits=NA, na.rm=TRUE, as.raster=FALSE) 

# S4 method for class 'SpatVector'
unique(x, incomparables=FALSE, geom=TRUE, atts=TRUE, ...)
```

## Arguments

- x:

  SpatRaster or SpatVector

- incomparables:

  logical. If `FALSE` and `x` is a SpatRaster: the unique values are
  determined for all layers together, and the result is a matrix. If
  `TRUE`, each layer is evaluated separately, and a list is returned. If
  `x` is a SpatVector this argument is as for a `data.frame`

- digits:

  integer. The number of digits for rounding the values before finding
  the unique values. Use `NA` means to not do any rounding

- na.rm:

  logical. If `TRUE`, `NaN` is included if there are any missing values

- as.raster:

  logical. If `TRUE`, a single-layer categorical SpatRaster with the
  unique values is returned

- ...:

  additional arguments passed on to
  [`unique`](https://rdrr.io/r/base/unique.html) or
  [`identical`](https://rdrr.io/r/base/identical.html)

- geom:

  logical. If `TRUE` the geometries are considered to determine
  uniqueness

- atts:

  logical. If `TRUE` the attribute values are considered to determine
  uniqueness

## Value

If `x` is a SpatRaster: data.frame or list (if `incomparables=FALSE`)

If `x` is a SpatVector: SpatVector

## Examples

``` r
r <- rast(ncols=5, nrows=5)
values(r) <- rep(1:5, each=5)
unique(r)
#>   lyr.1
#> 1     1
#> 2     2
#> 3     3
#> 4     4
#> 5     5
s <- c(r, round(r/3))
unique(s)
#>   lyr.1 lyr.1
#> 1     1     0
#> 2     2     1
#> 3     3     1
#> 4     4     1
#> 5     5     2
unique(s,TRUE)
#> [[1]]
#> [1] 1 2 3 4 5
#> 
#> [[2]]
#> [1] 0 1 2
#> 

unique(s, as.raster=TRUE)
#> class       : SpatRaster 
#> size        : 5, 5, 1  (nrow, ncol, nlyr)
#> resolution  : 72, 36  (x, y)
#> extent      : -180, 180, -90, 90  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (CRS84) (OGC:CRS84) 
#> source(s)   : memory
#> categories  : label, lyr.1, lyr.1.1 
#> name        : label 
#> min value   :   1_0 
#> max value   :   5_2 


v <- vect(cbind(x=c(1:5,1:5), y=c(5:1,5:1)), 
    crs="+proj=utm +zone=1 +datum=WGS84")
nrow(v)
#> [1] 10
u <- unique(v)
nrow(u)
#> [1] 5

values(v) <- c(1:5, 1:3, 5:4)
unique(v)
#>  class       : SpatVector 
#>  geometry    : points 
#>  dimensions  : 7, 1  (geometries, attributes)
#>  extent      : 1, 5, 1, 5  (xmin, xmax, ymin, ymax)
#>  coord. ref. : +proj=utm +zone=1 +datum=WGS84 +units=m +no_defs 
#>  names       : value
#>  type        : <int>
#>  values      :     1
#>                    2
#>                    3
```
