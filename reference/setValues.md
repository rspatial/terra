# Set the values of raster cells or of geometry attributes

Set cell values of a SpatRaster or the attributes of a SpatVector. For
large SpatRasters use
[`init`](https://rspatial.github.io/terra/reference/init.md) instead to
set values.

## Usage

``` r
# S4 method for class 'SpatRaster,ANY'
values(x) <- value

# S4 method for class 'SpatRaster,ANY'
setValues(x, values, keeptime=TRUE, keepunits=TRUE, keepnames=FALSE, props=FALSE)

# S4 method for class 'SpatVector,ANY'
values(x) <- value
```

## Arguments

- x:

  SpatRaster or SpatVector

- value:

  For SpatRaster: numeric, matrix or data.frame. The length of the
  numeric values must match the total number of cells (ncell(x) \*
  nlyr(x)), or be a single value. The number of columns of the matrix or
  data.frame must match the number of layers of `x`, and the number of
  rows must match the number of cells of `x`. It is also possible to use
  a matrix with the same number of rows as `x` and the number of columns
  that matches `ncol(x) * nlyr(x)`.

  For SpatVector: data.frame, matrix, vector, or NULL

- values:

  Same as for `value`

- keeptime:

  logical. If `TRUE` the time stamps are kept

- keepunits:

  logical. If `FALSE` the units are discarded

- keepnames:

  logical. If `FALSE` the layer names are replaced by the column names
  in `y` (if present)

- props:

  logical. If `TRUE` the properties (categories and color-table) are
  kept

## Value

The same object type as `x`

## See also

[`values`](https://rspatial.github.io/terra/reference/values.md),
[`init`](https://rspatial.github.io/terra/reference/init.md)

## Examples

``` r
f <- system.file("ex/elev.tif", package="terra")
r <- rast(f)
x <- setValues(r, 1:ncell(r))
x
#> class       : SpatRaster 
#> size        : 90, 95, 1  (nrow, ncol, nlyr)
#> resolution  : 0.008333333, 0.008333333  (x, y)
#> extent      : 5.741667, 6.533333, 49.44167, 50.19167  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (EPSG:4326) 
#> source(s)   : memory
#> varname     : elev 
#> name        : elevation 
#> min value   :         1 
#> max value   :      8550 
values(x) <- runif(ncell(x))
x
#> class       : SpatRaster 
#> size        : 90, 95, 1  (nrow, ncol, nlyr)
#> resolution  : 0.008333333, 0.008333333  (x, y)
#> extent      : 5.741667, 6.533333, 49.44167, 50.19167  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (EPSG:4326) 
#> source(s)   : memory
#> varname     : elev 
#> name        :    elevation 
#> min value   : 0.0001064336 
#> max value   : 0.9999305937 
head(x)
#>    elevation
#> 1 0.07706438
#> 2 0.03554058
#> 3 0.64279549
#> 4 0.92861520
#> 5 0.59809242
#> 6 0.56090075

f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
values(v) <- data.frame(ID=1:12, name=letters[1:12])
head(v)
#>   ID name
#> 1  1    a
#> 2  2    b
#> 3  3    c
#> 4  4    d
#> 5  5    e
#> 6  6    f
```
