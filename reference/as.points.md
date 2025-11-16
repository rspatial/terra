# Conversion to a SpatVector of points

Conversion of a SpatRaster, SpatVector or SpatExtent to a SpatVector of
points.

## Usage

``` r
# S4 method for class 'SpatRaster'
as.points(x, values=TRUE, na.rm=TRUE, na.all=FALSE)

# S4 method for class 'SpatVector'
as.points(x, multi=FALSE, skiplast=TRUE)

# S4 method for class 'SpatExtent'
as.points(x, crs="")
```

## Arguments

- x:

  SpatRaster, SpatVector or SpatExtent

- values:

  logical; include cell values as attributes?

- multi:

  logical. If `TRUE` a multi-point geometry is returned

- skiplast:

  logical. If `TRUE` the last point of a polygon (which is the same as
  the first point) is not included

- na.rm:

  logical. If `TRUE` cells that are `NA` are ignored

- na.all:

  logical. If `TRUE` cells are only ignored if `na.rm=TRUE` and their
  value is `NA` for **all** layers instead of for `any` layer

- crs:

  character. The coordinate reference system (see
  [`crs`](https://rspatial.github.io/terra/reference/crs.md))

## Value

SpatVector

## See also

[`as.lines`](https://rspatial.github.io/terra/reference/as.lines.md),
`as.points`

## Examples

``` r
r <- rast(ncols=2, nrows=2)
values(r) <- 1:ncell(r)

as.points(r)
#>  class       : SpatVector 
#>  geometry    : points 
#>  dimensions  : 4, 1  (geometries, attributes)
#>  extent      : -90, 90, -45, 45  (xmin, xmax, ymin, ymax)
#>  coord. ref. : lon/lat WGS 84 (CRS84) (OGC:CRS84) 
#>  names       : lyr.1
#>  type        : <num>
#>  values      :     1
#>                    2
#>                    3

p <- as.polygons(r)
as.points(p)
#>  class       : SpatVector 
#>  geometry    : points 
#>  dimensions  : 16, 1  (geometries, attributes)
#>  extent      : -180, 180, -90, 90  (xmin, xmax, ymin, ymax)
#>  coord. ref. : lon/lat WGS 84 (CRS84) (OGC:CRS84) 
#>  names       : lyr.1
#>  type        : <int>
#>  values      :     1
#>                    1
#>                    1
```
