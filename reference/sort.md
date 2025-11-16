# Sort a SpatRaster or SpatVector

Sort the cell values of a SpatRaster across layers. You can also compute
the sorting order.

Or sort the records of SpatVector (or data.frame) by specifying the
column number(s) or names(s) to sort on.

## Usage

``` r
# S4 method for class 'SpatRaster'
sort(x, decreasing=FALSE, order=FALSE, filename="", ...)

# S4 method for class 'SpatVector'
sort(x, v, decreasing=FALSE)
```

## Arguments

- x:

  SpatRaster

- decreasing:

  logical. If `TRUE`, sorting is in decreasing order

- order:

  logical. If `TRUE` the sorting order is returned instead of the sorted
  values

- filename:

  character. Output filename

- ...:

  additional arguments for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

- v:

  character or numeric indicating the column(s) to sort on

## Value

SpatRaster

## Examples

``` r
f <- system.file("ex/elev.tif", package="terra")
r <- rast(f)
r <- c(r, r/2, r*2)
sort(r)
#> class       : SpatRaster 
#> size        : 90, 95, 3  (nrow, ncol, nlyr)
#> resolution  : 0.008333333, 0.008333333  (x, y)
#> extent      : 5.741667, 6.533333, 49.44167, 50.19167  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (EPSG:4326) 
#> source(s)   : memory
#> varname     : elev 
#> names       : elevation, elevation, elevation 
#> min values  :      70.5,       141,       282 
#> max values  :     273.5,       547,      1094 

ord <- sort(r, order=TRUE)
# these two are the same 
ord[[1]]
#> class       : SpatRaster 
#> size        : 90, 95, 1  (nrow, ncol, nlyr)
#> resolution  : 0.008333333, 0.008333333  (x, y)
#> extent      : 5.741667, 6.533333, 49.44167, 50.19167  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (EPSG:4326) 
#> source(s)   : memory
#> varname     : elev 
#> name        : elevation 
#> min value   :         0 
#> max value   :         1 
which.min(r)
#> class       : SpatRaster 
#> size        : 90, 95, 1  (nrow, ncol, nlyr)
#> resolution  : 0.008333333, 0.008333333  (x, y)
#> extent      : 5.741667, 6.533333, 49.44167, 50.19167  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (EPSG:4326) 
#> source(s)   : memory
#> name        : which.min 
#> min value   :         2 
#> max value   :         2 
```
