# Create a SpatVectorCollection

Methods to create a SpatVectorCollection. This is an object to hold
"sub-datasets", each a SpatVector, perhaps of different geometry type.

## Usage

``` r
# S4 method for class 'missing'
svc(x) 

# S4 method for class 'SpatVector'
svc(x, ...) 

# S4 method for class 'list'
svc(x) 

# S4 method for class 'character'
svc(x, layer="", query="", dialect="", extent=NULL, filter=NULL)
```

## See also

[`sprc`](https://rspatial.github.io/terra/reference/sprc.md)

## Arguments

- x:

  SpatVector, character (filename), list with SpatVectors, or missing

- ...:

  Additional `SpatVector`s

- layer:

  character. layer name to select a layer from a file (database) with
  multiple layers

- query:

  character. A query to subset the dataset

- dialect:

  character. The SQL dialect to use (if any). For example: "SQLite". ""
  refers to the default [OGR-SQL
  dialect](https://gdal.org/en/latest/user/ogr_sql_dialect.html)

- extent:

  Spat\* object. The extent of the object is used as a spatial filter to
  select the geometries to read. Ignored if `filter` is not `NULL`

- filter:

  SpatVector. Used as a spatial filter to select geometries to read (the
  convex hull is used for lines or points). It is guaranteed that all
  features that overlap with the extent of filter will be returned. It
  can happen that additional geometries are returned

## Value

SpatVectorCollection

## Examples

``` r
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
x <- svc()
x <- svc(v, v[1:3,], as.lines(v[3:5,]), as.points(v))
length(x)
#> [1] 4
x
#>  class       : SpatVectorCollection 
#>  length      : 4 
#>  geometry    : polygons (12)
#>                polygons (3)
#>                lines (3)
#>                points (3983)
#>  crs (first) : lon/lat WGS 84 (EPSG:4326) 
#>  names       : , , ,  

# extract
x[3]
#>  class       : SpatVector 
#>  geometry    : lines 
#>  dimensions  : 3, 6  (geometries, attributes)
#>  extent      : 5.74414, 6.239243, 49.69933, 50.03632  (xmin, xmax, ymin, ymax)
#>  coord. ref. : lon/lat WGS 84 (EPSG:4326) 
#>  names       :  ID_1   NAME_1  ID_2  NAME_2  AREA       POP
#>  type        : <num>    <chr> <num>   <chr> <num>     <num>
#>  values      :     1 Diekirch     3 Redange   259 1.866e+04
#>                    1 Diekirch     4 Vianden    76      5163
#>                    1 Diekirch     5   Wiltz   263 1.674e+04

# replace
x[2] <- as.lines(v[1,])

```
