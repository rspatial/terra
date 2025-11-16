# Conversion to a SpatVector of polygons

Conversion of a SpatRaster, SpatVector or SpatExtent to a SpatVector of
polygons.

## Usage

``` r
# S4 method for class 'SpatRaster'
as.polygons(x, round=TRUE, aggregate=TRUE, values=TRUE,
    na.rm=TRUE, na.all=FALSE, extent=FALSE, digits=0, ...)

# S4 method for class 'SpatVector'
as.polygons(x, extent=FALSE)

# S4 method for class 'SpatExtent'
as.polygons(x, crs="")
```

## Arguments

- x:

  SpatRaster, SpatVector or SpatExtent

- round:

  logical; If `TRUE` and `aggregate=TRUE`, values are rounded before
  aggregation. If this value is `FALSE` the SpatVector returned can have
  very many polygons and can be very large

- aggregate:

  logical; combine cells with the same values? If `TRUE` only the first
  layer in `x` is processed

- values:

  logical; include cell values as attributes?

- extent:

  logical. if `TRUE`, a polygon for the extent of the SpatRaster or
  SpatVector is returned. If `x` is a SpatRaster, the polygon has
  vertices for each row and column, not just the four corners of the
  raster. This can be useful for more precise projection. If that is not
  required, it is more efficient to get the extent represented by only
  the four corners with `as.polygons(ext(x), crs=crs(x))`

- na.rm:

  logical. If `TRUE` cells that are `NA` are ignored

- na.all:

  logical. If `TRUE` cells are only ignored if `na.rm=TRUE` and their
  value is `NA` for **all** layers instead of for `any` layer

- digits:

  integer. The number of digits for rounding (if `round=TRUE`)

- crs:

  character. The coordinate reference system (see
  [`crs`](https://rspatial.github.io/terra/reference/crs.md))

- ...:

  additional arguments. For backward compatibility. Will be removed in
  the future

## Value

SpatVector

## See also

[`as.lines`](https://rspatial.github.io/terra/reference/as.lines.md),
[`as.points`](https://rspatial.github.io/terra/reference/as.points.md)

## Examples

``` r
r <- rast(ncols=2, nrows=2)
values(r) <- 1:ncell(r)

p <- as.polygons(r)
p
#>  class       : SpatVector 
#>  geometry    : polygons 
#>  dimensions  : 4, 1  (geometries, attributes)
#>  extent      : -180, 180, -90, 90  (xmin, xmax, ymin, ymax)
#>  coord. ref. : lon/lat WGS 84 (CRS84) (OGC:CRS84) 
#>  names       : lyr.1
#>  type        : <int>
#>  values      :     1
#>                    2
#>                    3
```
