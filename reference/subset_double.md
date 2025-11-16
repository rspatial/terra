# Subset a SpatRaster or a SpatVector

Select a subset of layers from a SpatRaster or select a subset of
records (row) and/or variables (columns) from a SpatVector.

## Usage

``` r
# S4 method for class 'SpatRaster,numeric,missing'
x[[i, j]]

# S4 method for class 'SpatRasterDataset,ANY,ANY'
x[[i, j, drop = TRUE]]

# S4 method for class 'SpatVector,numeric,missing'
x[[i, j, drop = FALSE]]
```

## Arguments

- x:

  SpatRaster or SpatVector

- i:

  if `x` is a `SpatRaster`: integer, logical, or character to select
  layers

  if `x` is a `SpatVector`: integer, logical, or character to select
  variables

- j:

  missing, or, for SpatRasterDataset only, numeric

- drop:

  logical. If `TRUE`, the geometries will be dropped, and a data.frame
  is returned

## Value

if `x` is a `SpatRaster` or `SpatRasterDataset`: SpatRaster

if `x` is a `SpatVector`: a `data.frame`.

## See also

[`subset`](https://rspatial.github.io/terra/reference/subset.md),
[`$`](https://rspatial.github.io/terra/reference/subset_dollar.md),
[`[`](https://rspatial.github.io/terra/reference/subset_single.md),
[`extract`](https://rspatial.github.io/terra/reference/extract.md)

## Examples

``` r
### SpatRaster
s <- rast(system.file("ex/logo.tif", package="terra"))   
s[[ 1:2 ]]
#> class       : SpatRaster 
#> size        : 77, 101, 2  (nrow, ncol, nlyr)
#> resolution  : 1, 1  (x, y)
#> extent      : 0, 101, 0, 77  (xmin, xmax, ymin, ymax)
#> coord. ref. : Cartesian (Meter) 
#> source      : logo.tif 
#> names       : red, green 
#> min values  :   0,     0 
#> max values  : 255,   255 

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

# does not with double brackets
# s[["re"]]


### SpatVector

v <- vect(system.file("ex/lux.shp", package="terra"))
v[[2:3]]
#>          NAME_1 ID_2
#> 1      Diekirch    1
#> 2      Diekirch    2
#> 3      Diekirch    3
#> 4      Diekirch    4
#> 5      Diekirch    5
#> 6  Grevenmacher    6
#> 7  Grevenmacher    7
#> 8  Grevenmacher   12
#> 9    Luxembourg    8
#> 10   Luxembourg    9
#> 11   Luxembourg   10
#> 12   Luxembourg   11

# to keep the geometry use
v[,2:3]
#>  class       : SpatVector 
#>  geometry    : polygons 
#>  dimensions  : 12, 2  (geometries, attributes)
#>  extent      : 5.74414, 6.528252, 49.44781, 50.18162  (xmin, xmax, ymin, ymax)
#>  source      : lux.shp
#>  coord. ref. : lon/lat WGS 84 (EPSG:4326) 
#>  names       :   NAME_1  ID_2
#>  type        :    <chr> <num>
#>  values      : Diekirch     1
#>                Diekirch     2
#>                Diekirch     3
```
