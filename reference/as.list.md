# Coerce a Spat\* object to a list

Coerce a SpatRaster, SpatRasterCollection, SpatRasterDataset, SpatVector
or SpatVectorCollection to a list. With a SpatRaster, each layer becomes
a list element. With a SpatRasterCollection or SpatRasterDataset, each
SpatRaster becomes a list element. With a SpatVector, each variable
(attribute) becomes a list element. With a SpatVectorCollection, each
SpatVector becomes a list element.

## Usage

``` r
# S4 method for class 'SpatRaster'
as.list(x, geom=NULL, ...)

# S4 method for class 'SpatRasterCollection'
as.list(x, ...)

# S4 method for class 'SpatVector'
as.list(x, geom=NULL, ...)

# S4 method for class 'SpatVectorCollection'
as.list(x, ...)
```

## Arguments

- x:

  SpatRaster, SpatRasterDataset, SpatRasterCollection, or SpatVector

- geom:

  character or NULL. If not NULL, and `x` is a SpatVector, it should be
  either "WKT" or "HEX", to get the geometry included in Well-Known-Text
  or hexadecimal notation. If `x` has point geometry, it can also bey
  "XY" to add the coordinates of each point. If `x` is a SpatRaster, any
  value that is not NULL will return a list with the the parameters
  describing the geometry of the SpatRaster are returned

- ...:

  additional arguments. These are ignored

## See also

see `coerce` for `as.data.frame` with a SpatRaster; and
[`geom`](https://rspatial.github.io/terra/reference/geometry.md) to only
extract the geometry of a SpatVector

## Value

list

## Examples

``` r
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
as.list(v)
#> $ID_1
#>  [1] 1 1 1 1 1 2 2 2 3 3 3 3
#> 
#> $NAME_1
#>  [1] "Diekirch"     "Diekirch"     "Diekirch"     "Diekirch"     "Diekirch"    
#>  [6] "Grevenmacher" "Grevenmacher" "Grevenmacher" "Luxembourg"   "Luxembourg"  
#> [11] "Luxembourg"   "Luxembourg"  
#> 
#> $ID_2
#>  [1]  1  2  3  4  5  6  7 12  8  9 10 11
#> 
#> $NAME_2
#>  [1] "Clervaux"         "Diekirch"         "Redange"          "Vianden"         
#>  [5] "Wiltz"            "Echternach"       "Remich"           "Grevenmacher"    
#>  [9] "Capellen"         "Esch-sur-Alzette" "Luxembourg"       "Mersch"          
#> 
#> $AREA
#>  [1] 312 218 259  76 263 188 129 210 185 251 237 233
#> 
#> $POP
#>  [1]  18081  32543  18664   5163  16735  18899  22366  29828  48187 176820
#> [11] 182607  32112
#> 


s <- rast(system.file("ex/logo.tif", package="terra")) + 1  
as.list(s)
#> [[1]]
#> class       : SpatRaster 
#> size        : 77, 101, 1  (nrow, ncol, nlyr)
#> resolution  : 1, 1  (x, y)
#> extent      : 0, 101, 0, 77  (xmin, xmax, ymin, ymax)
#> coord. ref. : Cartesian (Meter) 
#> source(s)   : memory
#> varname     : logo 
#> name        : red 
#> min value   :   1 
#> max value   : 256 
#> 
#> [[2]]
#> class       : SpatRaster 
#> size        : 77, 101, 1  (nrow, ncol, nlyr)
#> resolution  : 1, 1  (x, y)
#> extent      : 0, 101, 0, 77  (xmin, xmax, ymin, ymax)
#> coord. ref. : Cartesian (Meter) 
#> source(s)   : memory
#> varname     : logo 
#> name        : green 
#> min value   :     1 
#> max value   :   256 
#> 
#> [[3]]
#> class       : SpatRaster 
#> size        : 77, 101, 1  (nrow, ncol, nlyr)
#> resolution  : 1, 1  (x, y)
#> extent      : 0, 101, 0, 77  (xmin, xmax, ymin, ymax)
#> coord. ref. : Cartesian (Meter) 
#> source(s)   : memory
#> varname     : logo 
#> name        : blue 
#> min value   :    1 
#> max value   :  256 
#> 
```
