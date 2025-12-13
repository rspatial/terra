# SpatVector geometric properties

`width` returns the minimum diameter of the geometry, defined as the
smallest band that contains the geometry, where a band is a strip of the
plane defined by two parallel lines. This can be thought of as the
smallest hole that the geometry can be moved through, with a single
rotation.

`clearance` returns the minimum clearance of a geometry. The minimum
clearance is the smallest amount by which a vertex could be moved to
produce an invalid polygon, a non-simple linestring, or a multipoint
with repeated points. If a geometry has a minimum clearance of 'mc', it
can be said that "no two distinct vertices in the geometry are separated
by less than "mc". No vertex is closer than "mc" to a line segment of
which it is not an endpoint".

If the minimum clearance cannot be defined for a geometry (such as with
a single point), `NA` is returned.

## Usage

``` r
# S4 method for class 'SpatVector'
width(x, as.lines=FALSE)
# S4 method for class 'SpatVector'
clearance(x, as.lines=FALSE)
```

## Arguments

- x:

  SpatVector of lines or polygons

- as.lines:

  logical. If `TRUE` lines are returned that define the width or
  clearance

## Value

numeric or SpatVector

## See also

[`hull`](https://rspatial.github.io/terra/reference/convhull.md)

## Examples

``` r
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)

width(v)
#>  [1] 25511.501 20333.498 19019.103  9255.223 20141.189 14308.518 15475.916
#>  [8] 18111.001 17435.502 15750.096 18936.122 15692.522
clearance(v)
#>  [1] 24.349993 50.070894 31.878312 51.390103 13.517168 15.960117  2.871879
#>  [8] 17.584083  9.023613  4.239727  8.231322 60.918397

width_lines <- width(v, as.lines=TRUE)
plot(v)
lines(width_lines, col="blue")
```
