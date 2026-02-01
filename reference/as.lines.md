# Conversion to a SpatVector of lines

Conversion of a SpatRaster, SpatVector or SpatExtent to a SpatVector of
lines.

## Usage

``` r
# S4 method for class 'SpatRaster'
as.lines(x)

# S4 method for class 'SpatVector'
as.lines(x)

# S4 method for class 'SpatExtent'
as.lines(x, crs="")

# S4 method for class 'matrix'
as.lines(x, crs="", segments=FALSE)
```

## Arguments

- x:

  SpatRaster, SpatVector, SpatExtent or matrix. If `x` is a matrix it
  should have two columns for a single line, or four columns, where each
  row has the start and end coordinates (x, y) for lines

- crs:

  character. The coordinate reference system (see
  [`crs`](https://rspatial.github.io/terra/reference/crs.md))

- segments:

  logical. Should (poly-)lines or polygons be disaggregated into their
  line-segments? See
  [`disagg`](https://rspatial.github.io/terra/reference/disaggregate.md)

## Value

SpatVector

## See also

[`as.points`](https://rspatial.github.io/terra/reference/as.points.md),
[`as.polygons`](https://rspatial.github.io/terra/reference/as.polygons.md)

## Examples

``` r
r <- rast(ncols=2, nrows=2)
values(r) <- 1:ncell(r)

as.lines(r)
#>  class       : SpatVector 
#>  geometry    : lines 
#>  dimensions  : 6, 0  (geometries, attributes)
#>  extent      : -180, 180, -90, 90  (xmin, xmax, ymin, ymax)
#>  coord. ref. : lon/lat WGS 84 (CRS84) (OGC:CRS84) 

as.lines(ext(r), crs=crs(r))
#>  class       : SpatVector 
#>  geometry    : lines 
#>  dimensions  : 1, 0  (geometries, attributes)
#>  extent      : -180, 180, -90, 90  (xmin, xmax, ymin, ymax)
#>  coord. ref. : lon/lat WGS 84 (CRS84) (OGC:CRS84) 

p <- as.polygons(r)
as.lines(p)
#>  class       : SpatVector 
#>  geometry    : lines 
#>  dimensions  : 4, 1  (geometries, attributes)
#>  extent      : -180, 180, -90, 90  (xmin, xmax, ymin, ymax)
#>  coord. ref. : lon/lat WGS 84 (CRS84) (OGC:CRS84) 
#>  names       : lyr.1
#>  type        : <int>
#>  values      :     1
#>                    2
#>                    3


## with a matrix
s <- cbind(1:5, 1:5)
e <- cbind(1:5, 0)

as.lines(s)
#>  class       : SpatVector 
#>  geometry    : lines 
#>  dimensions  : 1, 0  (geometries, attributes)
#>  extent      : 1, 5, 1, 5  (xmin, xmax, ymin, ymax)
#>  coord. ref. :  
as.lines(cbind(s, e), "+proj=longlat")
#>  class       : SpatVector 
#>  geometry    : lines 
#>  dimensions  : 5, 0  (geometries, attributes)
#>  extent      : 1, 5, 0, 5  (xmin, xmax, ymin, ymax)
#>  coord. ref. : +proj=longlat +datum=WGS84 +no_defs 
```
