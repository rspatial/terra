# Extract values from a SpatRaster, SpatVector or SpatExtent

Extract values from a SpatRaster; a subset of records (row) and/or
variables (columns) from a SpatVector; or a number from a SpatExtent.

You can use indices (row, column, layer or cell numbers) to extract. You
can also use other Spat\* objects.

## Usage

``` r
# S4 method for class 'SpatRaster,ANY,ANY,ANY'
x[i, j, k]

# S4 method for class 'SpatVector,numeric,numeric'
x[i, j, drop = FALSE]

# S4 method for class 'SpatVector,SpatVector,missing'
x[i, j]

# S4 method for class 'SpatExtent,numeric,missing'
x[i, j]
```

## Arguments

- x:

  SpatRaster, SpatVector or SpatExtent

- i:

  if `x` is a `SpatRaster`: numeric, logical or missing to select rows
  or, if `j` is missing, to select cells numbers.

  if `x` is a `SpatVector`: numeric or missing to select rows. if `i` is
  another SpatVector: get a new SpatVector with the geometries that
  intersect.

  if `x` is a `SpatExtent`: integer between 1 and 4.

- j:

  numeric, logical, or missing to select columns

- k:

  numeric, character, or missing to select layers

- drop:

  logical. If `FALSE` an object of the same class as `x` is returned

## Value

numeric if `x` is a SpatExtent. Same as `x` if `drop=FALSE`. Otherwise a
data.frame

## See also

[`extract`](https://rspatial.github.io/terra/reference/extract.md),
[`subset`](https://rspatial.github.io/terra/reference/subset.md),
[`$`](https://rspatial.github.io/terra/reference/subset_dollar.md),
[`[[`](https://rspatial.github.io/terra/reference/subset_double.md)

## Examples

``` r
### SpatRaster
f <- system.file("ex/elev.tif", package="terra")
r <- rast(f)
r[3638]
#>   elevation
#> 1       473
rowColFromCell(r, 2638)
#>      [,1] [,2]
#> [1,]   28   73
r[39, 28]
#>   elevation
#> 1       473
x <- r[39:40, 28:29, drop=FALSE] 
as.matrix(x, wide=TRUE)
#>      [,1] [,2]
#> [1,]  473  429
#> [2,]  430  436

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
v[1:2, 2:3, drop=TRUE]
#>     NAME_1 ID_2
#> 1 Diekirch    1
#> 2 Diekirch    2
```
