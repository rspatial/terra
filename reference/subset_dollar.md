# Subset a SpatRaster or a SpatVector

Select a subset of layers from a SpatRaster or select a subset of
records (row) and/or variables (columns) from a SpatVector.

## Usage

``` r
# S4 method for class 'SpatExtent'
x$name
```

## Arguments

- x:

  SpatRaster, SpatVector or SpatExtent

- name:

  character. If `x` is a SpatRaster: layer name. If `x` is a SpatVector:
  variable name. If `x` is a SpatExtent: xmin, xmax, ymin or ymax

## Value

if `x` is a `SpatRaster`: SpatRaster

if `x` is a `SpatVector`: SpatVector or, if `drop=TRUE`, a `data.frame`.

## See also

[`subset`](https://rspatial.github.io/terra/reference/subset.md),
[`[`](https://rspatial.github.io/terra/reference/subset_single.md),
[`[[`](https://rspatial.github.io/terra/reference/subset_double.md),
[`extract`](https://rspatial.github.io/terra/reference/extract.md)

## Examples

``` r
### SpatRaster
s <- rast(system.file("ex/logo.tif", package="terra"))   
subset(s, 2:3)
#> class       : SpatRaster 
#> size        : 77, 101, 2  (nrow, ncol, nlyr)
#> resolution  : 1, 1  (x, y)
#> extent      : 0, 101, 0, 77  (xmin, xmax, ymin, ymax)
#> coord. ref. : Cartesian (Meter) 
#> source      : logo.tif 
#> names       : green, blue 
#> min values  :     0,    0 
#> max values  :   255,  255 
subset(s, c(3,2,3,1))
#> class       : SpatRaster 
#> size        : 77, 101, 4  (nrow, ncol, nlyr)
#> resolution  : 1, 1  (x, y)
#> extent      : 0, 101, 0, 77  (xmin, xmax, ymin, ymax)
#> coord. ref. : Cartesian (Meter) 
#> source      : logo.tif 
#> names       : blue, green, blue, red 
#> min values  :    0,     0,    0,   0 
#> max values  :  255,   255,  255, 255 
#equivalent to 
s[[ c(3,2,3,1) ]]
#> class       : SpatRaster 
#> size        : 77, 101, 4  (nrow, ncol, nlyr)
#> resolution  : 1, 1  (x, y)
#> extent      : 0, 101, 0, 77  (xmin, xmax, ymin, ymax)
#> coord. ref. : Cartesian (Meter) 
#> source      : logo.tif 
#> names       : blue, green, blue, red 
#> min values  :    0,     0,    0,   0 
#> max values  :  255,   255,  255, 255 

s[[c("red", "green")]]
#> class       : SpatRaster 
#> size        : 77, 101, 2  (nrow, ncol, nlyr)
#> resolution  : 1, 1  (x, y)
#> extent      : 0, 101, 0, 77  (xmin, xmax, ymin, ymax)
#> coord. ref. : Cartesian (Meter) 
#> source      : logo.tif 
#> names       : red, green 
#> min values  :   0,     0 
#> max values  : 255,   255 
s$red
#> class       : SpatRaster 
#> size        : 77, 101, 1  (nrow, ncol, nlyr)
#> resolution  : 1, 1  (x, y)
#> extent      : 0, 101, 0, 77  (xmin, xmax, ymin, ymax)
#> coord. ref. : Cartesian (Meter) 
#> source      : logo.tif 
#> name        : red 
#> min value   :   0 
#> max value   : 255 

# expression based (partial) matching of names with single brackets
s["re"]
#> class       : SpatRaster 
#> size        : 77, 101, 2  (nrow, ncol, nlyr)
#> resolution  : 1, 1  (x, y)
#> extent      : 0, 101, 0, 77  (xmin, xmax, ymin, ymax)
#> coord. ref. : Cartesian (Meter) 
#> source      : logo.tif 
#> names       : red, green 
#> min values  :   0,     0 
#> max values  : 255,   255 
s["^re"]
#> class       : SpatRaster 
#> size        : 77, 101, 1  (nrow, ncol, nlyr)
#> resolution  : 1, 1  (x, y)
#> extent      : 0, 101, 0, 77  (xmin, xmax, ymin, ymax)
#> coord. ref. : Cartesian (Meter) 
#> source      : logo.tif 
#> name        : red 
#> min value   :   0 
#> max value   : 255 

# not with double brackets
# s[["re"]]


### SpatVector

v <- vect(system.file("ex/lux.shp", package="terra"))
v[2:3,]
#>  class       : SpatVector 
#>  geometry    : polygons 
#>  dimensions  : 2, 6  (geometries, attributes)
#>  extent      : 5.746118, 6.315773, 49.69933, 49.98745  (xmin, xmax, ymin, ymax)
#>  coord. ref. : lon/lat WGS 84 (EPSG:4326) 
#>  names       :  ID_1   NAME_1  ID_2   NAME_2  AREA       POP
#>  type        : <num>    <chr> <num>    <chr> <num>     <num>
#>  values      :     1 Diekirch     2 Diekirch   218 3.254e+04
#>                    1 Diekirch     3  Redange   259 1.866e+04
v[1:2, 2:3]
#>  class       : SpatVector 
#>  geometry    : polygons 
#>  dimensions  : 2, 2  (geometries, attributes)
#>  extent      : 5.826232, 6.315773, 49.78479, 50.18162  (xmin, xmax, ymin, ymax)
#>  coord. ref. : lon/lat WGS 84 (EPSG:4326) 
#>  names       :   NAME_1  ID_2
#>  type        :    <chr> <num>
#>  values      : Diekirch     1
#>                Diekirch     2

subset(v, v$NAME_1 == "Diekirch", c("NAME_1", "NAME_2"))
#>  class       : SpatVector 
#>  geometry    : polygons 
#>  dimensions  : 5, 2  (geometries, attributes)
#>  extent      : 5.74414, 6.315773, 49.69933, 50.18162  (xmin, xmax, ymin, ymax)
#>  coord. ref. : lon/lat WGS 84 (EPSG:4326) 
#>  names       :   NAME_1   NAME_2
#>  type        :    <chr>    <chr>
#>  values      : Diekirch Clervaux
#>                Diekirch Diekirch
#>                Diekirch  Redange

subset(v, NAME_1 == "Diekirch", c(NAME_1, NAME_2), NSE=TRUE)
#>  class       : SpatVector 
#>  geometry    : polygons 
#>  dimensions  : 5, 2  (geometries, attributes)
#>  extent      : 5.74414, 6.315773, 49.69933, 50.18162  (xmin, xmax, ymin, ymax)
#>  coord. ref. : lon/lat WGS 84 (EPSG:4326) 
#>  names       :   NAME_1   NAME_2
#>  type        :    <chr>    <chr>
#>  values      : Diekirch Clervaux
#>                Diekirch Diekirch
#>                Diekirch  Redange
```
