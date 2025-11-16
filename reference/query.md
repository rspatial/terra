# Query a SpatVectorProxy object

Query a SpatVectorProxy to extract a subset

## Usage

``` r
# S4 method for class 'SpatVectorProxy'
query(x, start=1, n=nrow(x), vars=NULL, where=NULL, 
    extent=NULL, filter=NULL, sql=NULL, dialect="", what="")
```

## Arguments

- x:

  SpatVectorProxy

- start:

  positive integer. The record to start reading at

- n:

  positive integer. The number of records requested

- vars:

  character. Variable names. Must be a subset of `names(x)`

- where:

  character. expression like "NAME_1='California' AND ID \> 3" , to
  subset records. Note that start and n are applied after executing the
  where statement

- extent:

  Spat\* object. The extent of the object is used as a spatial filter to
  select the geometries to read. Ignored if `filter` is not `NULL`

- filter:

  SpatVector. Used as a spatial filter to select geometries to read (the
  convex hull is used for lines or points)

- sql:

  character. Arbitrary SQL statement. If used, arguments "start", "n",
  "vars" and "where" are ignored

- what:

  character indicating what to read. Either `""` for geometries and
  attributes, or `"geoms"` to only read the geometries, `"attributes"`
  to only read the attributes (that are returned as a data.frame)

- dialect:

  character. The SQL dialect to use (if any). For example: "SQLite"

## See also

[`vect`](https://rspatial.github.io/terra/reference/vect.md)

## Value

SpatVector

## Examples

``` r
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f, proxy=TRUE)
v
#>  class       : SpatVectorProxy
#>  geometry    : polygons 
#>  dimensions  : 12, 6  (geometries, attributes)
#>  extent      : 5.74414, 6.528252, 49.44781, 50.18162  (xmin, xmax, ymin, ymax)
#>  source      : lux.shp
#>  coord. ref. : lon/lat WGS 84 (EPSG:4326) 
#>  names       :  ID_1 NAME_1  ID_2 NAME_2  AREA   POP
#>  type        : <num>  <chr> <num>  <chr> <num> <num>

x <- query(v, vars=c("ID_2", "NAME_2"), start=5, n=2)
x
#>  class       : SpatVector 
#>  geometry    : polygons 
#>  dimensions  : 2, 2  (geometries, attributes)
#>  extent      : 5.74414, 6.528252, 49.72324, 50.03632  (xmin, xmax, ymin, ymax)
#>  source      : lux.shp
#>  coord. ref. : lon/lat WGS 84 (EPSG:4326) 
#>  names       :  ID_2     NAME_2
#>  type        : <num>      <chr>
#>  values      :     5      Wiltz
#>                    6 Echternach

query(v, vars=c("ID_2", "NAME_1", "NAME_2"), where="NAME_1='Grevenmacher' AND ID_2 > 6")
#>  class       : SpatVector 
#>  geometry    : polygons 
#>  dimensions  : 2, 3  (geometries, attributes)
#>  extent      : 6.169137, 6.516485, 49.46498, 49.75016  (xmin, xmax, ymin, ymax)
#>  source      : lux.shp
#>  coord. ref. : lon/lat WGS 84 (EPSG:4326) 
#>  names       :  ID_2       NAME_1       NAME_2
#>  type        : <num>        <chr>        <chr>
#>  values      :     7 Grevenmacher       Remich
#>                   12 Grevenmacher Grevenmacher

## with an extent
e <- ext(5.9, 6.3, 49.9, 50) 
x <- query(v, extent=e)

## with polygons
p <- as.polygons(e)
x <- query(v, filter=p)
x
#>  class       : SpatVector 
#>  geometry    : polygons 
#>  dimensions  : 4, 6  (geometries, attributes)
#>  extent      : 5.74414, 6.315773, 49.78479, 50.18162  (xmin, xmax, ymin, ymax)
#>  source      : lux.shp
#>  coord. ref. : lon/lat WGS 84 (EPSG:4326) 
#>  names       :  ID_1   NAME_1  ID_2   NAME_2  AREA       POP
#>  type        : <num>    <chr> <num>    <chr> <num>     <num>
#>  values      :     1 Diekirch     1 Clervaux   312 1.808e+04
#>                    1 Diekirch     2 Diekirch   218 3.254e+04
#>                    1 Diekirch     4  Vianden    76      5163
```
