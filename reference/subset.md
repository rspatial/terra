# Subset a SpatRaster or a SpatVector

Select a subset of layers from a SpatRaster or select a subset of
records (row) and/or variables (columns) from a SpatVector.

## Usage

``` r
# S4 method for class 'SpatRaster'
subset(x, subset, negate=FALSE, NSE=FALSE, filename="", overwrite=FALSE, ...)

# S4 method for class 'SpatVector'
subset(x, subset, select, drop=FALSE, NSE=FALSE)
```

## Arguments

- x:

  SpatRaster or SpatVector

- subset:

  if `x` is a `SpatRaster`: integer or character to select layers

  if `x` is a `SpatVector`: logical expression indicating the rows to
  keep (missing values are taken as FALSE), or another `Spat*` object in
  which case the extent is used to spatially subset the intersecting
  geometries

- select:

  expression, indicating columns to select

- negate:

  logical. If `TRUE` all layers that are **not** in the subset are
  selected

- NSE:

  logical. If `TRUE`, non-standard evaluation (the use of unquoted
  variable names) is allowed. Set this to `FALSE` when calling `subset`
  from a function

- drop:

  logical. If `TRUE`, the geometries will be dropped, and a data.frame
  is returned

- filename:

  character. Output filename

- overwrite:

  logical. If `TRUE`, `filename` is overwritten

- ...:

  additional arguments for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## Value

if `x` is a `SpatRaster`: SpatRaster

if `x` is a `SpatVector`: SpatVector or, if `drop=TRUE`, a `data.frame`.

## See also

[`$`](https://rspatial.github.io/terra/reference/subset_dollar.md),
[`[[`](https://rspatial.github.io/terra/reference/subset_double.md),
[`[`](https://rspatial.github.io/terra/reference/subset_single.md)

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

# or like this
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
v[1:2, c("NAME_1", "NAME_2")]
#>  class       : SpatVector 
#>  geometry    : polygons 
#>  dimensions  : 2, 2  (geometries, attributes)
#>  extent      : 5.826232, 6.315773, 49.78479, 50.18162  (xmin, xmax, ymin, ymax)
#>  coord. ref. : lon/lat WGS 84 (EPSG:4326) 
#>  names       :   NAME_1   NAME_2
#>  type        :    <chr>    <chr>
#>  values      : Diekirch Clervaux
#>                Diekirch Diekirch

# or by location, i.e. by spatial overlap with another object
poly <- as.polygons(ext(6, 6.4, 49.75, 50))
subset(v, poly)
#>  class       : SpatVector 
#>  geometry    : polygons 
#>  dimensions  : 8, 6  (geometries, attributes)
#>  extent      : 5.74414, 6.528252, 49.58699, 50.18162  (xmin, xmax, ymin, ymax)
#>  coord. ref. : lon/lat WGS 84 (EPSG:4326) 
#>  names       :  ID_1   NAME_1  ID_2   NAME_2  AREA       POP
#>  type        : <num>    <chr> <num>    <chr> <num>     <num>
#>  values      :     1 Diekirch     1 Clervaux   312 1.808e+04
#>                    1 Diekirch     2 Diekirch   218 3.254e+04
#>                    1 Diekirch     3  Redange   259 1.866e+04
v[poly]
#>  class       : SpatVector 
#>  geometry    : polygons 
#>  dimensions  : 8, 6  (geometries, attributes)
#>  extent      : 5.74414, 6.528252, 49.58699, 50.18162  (xmin, xmax, ymin, ymax)
#>  coord. ref. : lon/lat WGS 84 (EPSG:4326) 
#>  names       :  ID_1   NAME_1  ID_2   NAME_2  AREA       POP
#>  type        : <num>    <chr> <num>    <chr> <num>     <num>
#>  values      :     1 Diekirch     1 Clervaux   312 1.808e+04
#>                    1 Diekirch     2 Diekirch   218 3.254e+04
#>                    1 Diekirch     3  Redange   259 1.866e+04
```
