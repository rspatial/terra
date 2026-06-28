# Distance to the furthest location on a geometry

For each point in `x`, compute the distance to the *furthest* location
on the edge of the geometries in `y`. This is the counterpart of
[`distance`](https://rspatial.github.io/terra/reference/distance.md),
which returns the distance to the *nearest* location.

The "edge" is the boundary of polygons, the line itself for lines, and
the points for a points geometry.

## Usage

``` r
# S4 method for class 'SpatVector,SpatVector'
furdist(x, y, pairwise=FALSE, unit="m")
```

## Arguments

- x:

  SpatVector of points. The location(s) from which the distance is
  measured

- y:

  SpatVector of points, lines, or polygons. The geometries whose
  furthest edge location is measured

- pairwise:

  logical. If `TRUE` and if `x` and `y` have the same number of
  geometries (or one of them has a single geometry), the pairwise
  distances (point `i` of `x` with geometry `i` of `y`) are returned. If
  `FALSE`, the distances for all combinations are returned

- unit:

  character. Can be either "m" or "km"

## Value

If `pairwise=TRUE` a numeric vector with one distance for each pair.
Otherwise a data.frame with columns "from" (the index of the point in
`x`), "to" (the index of the geometry in `y`) and "distance".

## See also

[`distance`](https://rspatial.github.io/terra/reference/distance.md),
[`nearest`](https://rspatial.github.io/terra/reference/nearby.md),
[`snapTo`](https://rspatial.github.io/terra/reference/snapTo.md)

## Examples

``` r
# planar
crs <- "+proj=utm +zone=1"
p <- vect("POLYGON ((0 0, 8 0, 8 9, 0 9, 0 0))", crs=crs)
pts <- vect(cbind(c(1, 5), c(1, 5)), crs=crs)

furdist(pts, p)
#>   from to  distance
#> 1    1  1 10.630146
#> 2    2  1  7.071068
furdist(pts, p, unit="km")
#>   from to    distance
#> 1    1  1 0.010630146
#> 2    2  1 0.007071068

# longitude/latitude
pp <- vect("POLYGON ((0 0, 10 0, 10 10, 0 10, 0 0))", crs="lonlat")
q <- vect(cbind(5, 5), crs="lonlat")
furdist(q, pp, unit="km")
#>   from to distance
#> 1    1  1 784.0287
```
