# SpatRaster or SpatVector to data.frame

Coerce a SpatRaster or SpatVector to a data.frame

## Usage

``` r
# S4 method for class 'SpatVector'
as.data.frame(x, row.names=NULL, optional=FALSE, geom=NULL, ...)

# S4 method for class 'SpatRaster'
as.data.frame(x, row.names=NULL, optional=FALSE, xy=FALSE, 
  cells=FALSE, time=FALSE, na.rm=NA, wide=TRUE, ...)
```

## Arguments

- x:

  SpatRaster or SpatVector

- geom:

  character or NULL. If not NULL, either "WKT" or "HEX", to get the
  geometry included in Well-Known-Text or hexadecimal notation. If `x`
  has point geometry, it can also be "XY" to add the coordinates of each
  point

- xy:

  logical. If `TRUE`, the coordinates of each raster cell are included

- time:

  logical. If `TRUE`, the time data is included (if available)

- na.rm:

  logical. If `TRUE`, cells that have a `NA` value in at least one layer
  are removed. If the argument is set to `NA` only cells that have `NA`
  values in all layers are removed

- cells:

  logical. If `TRUE`, the cell numbers of each raster cell are included

- wide:

  logical. If `FALSE`, the data.frame returned has a "long" format

- ...:

  Additional arguments passed to the
  [`data.frame`](https://rdrr.io/r/base/data.frame.html)

- row.names:

  This argument is ignored

- optional:

  This argument is ignored

## See also

[`as.list`](https://rspatial.github.io/terra/reference/as.list.md)`, `[`as.matrix`](https://rspatial.github.io/terra/reference/coerce.md).
See [`geom`](https://rspatial.github.io/terra/reference/geometry.md) to
only extract the geometry of a SpatVector

## Value

data.frame

## Examples

``` r
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
as.data.frame(v)
#>    ID_1       NAME_1 ID_2           NAME_2 AREA    POP
#> 1     1     Diekirch    1         Clervaux  312  18081
#> 2     1     Diekirch    2         Diekirch  218  32543
#> 3     1     Diekirch    3          Redange  259  18664
#> 4     1     Diekirch    4          Vianden   76   5163
#> 5     1     Diekirch    5            Wiltz  263  16735
#> 6     2 Grevenmacher    6       Echternach  188  18899
#> 7     2 Grevenmacher    7           Remich  129  22366
#> 8     2 Grevenmacher   12     Grevenmacher  210  29828
#> 9     3   Luxembourg    8         Capellen  185  48187
#> 10    3   Luxembourg    9 Esch-sur-Alzette  251 176820
#> 11    3   Luxembourg   10       Luxembourg  237 182607
#> 12    3   Luxembourg   11           Mersch  233  32112
```
