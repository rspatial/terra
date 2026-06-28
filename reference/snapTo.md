# Snap points to the edge of lines or polygons

Move the points in `x` onto the edge of the corresponding line or
polygon in `y`. A point that touches the the target geometry is not
changed. A point that is outside the geometry is moved to the nearest
location on the edge (the boundary of a polygon, or the line itself).

## Usage

``` r
# S4 method for class 'SpatVector,SpatVector'
snapTo(x, y, paired=TRUE, ...)
```

## Arguments

- x:

  SpatVector of points

- y:

  SpatVector of lines or polygons. Either with the same number of
  geometries as there are points in `x`, or with a single geometry

- paired:

  logical. If `TRUE` the points in `x` are paired with the geometries in
  `y` (point `i` with geometry `i`). `y` may instead have a single
  geometry, in which case all points are snapped to that geometry.
  Otherwise the points in `x` are snapped to any geometry in `y`

- ...:

  additional arguments. None implemented

## Value

SpatVector of points

## See also

[`nearest`](https://rspatial.github.io/terra/reference/nearby.md),
[`distance`](https://rspatial.github.io/terra/reference/distance.md),
[`furdist`](https://rspatial.github.io/terra/reference/furdist.md)

## Examples

``` r
crs <- "+proj=utm +zone=1"
p <- vect("POLYGON ((0 0, 8 0, 8 9, 0 9, 0 0))", crs=crs)

# one point inside the polygon (unchanged), one outside (snapped to the edge)
pts <- vect(cbind(c(4, 12), c(4, 4)), crs=crs)
s <- snapTo(pts, p)
crds(s)
#>      x y
#> [1,] 4 4
#> [2,] 8 4

# snap to a line
L <- vect("LINESTRING (0 0, 10 0)", crs=crs)
q <- vect(cbind(c(5, 5), c(0, 3)), crs=crs)
crds(snapTo(q, L))
#>      x y
#> [1,] 5 0
#> [2,] 5 0
```
