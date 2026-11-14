# Convex, concave, rectangular and circular hulls

Compute hulls around SpatVector geometries. This can be the convex hull,
the minimal bounding rotated rectangle, the minimal bounding circle, or
a concave hull. The concaveness of the concave hull can be specified in
different ways.

The old method `convHull` is deprecated and will be removed in a future
version.

## Usage

``` r
# S4 method for class 'SpatVector'
hull(x, type="convex", by="", param=1, allowHoles=TRUE, tight=TRUE, percent=100)
```

## Arguments

- x:

  SpatVector

- type:

  character. One of "convex", "rectangle", "circle", "concave_ratio",
  "concave_length"

- by:

  character (variable name), to get a new geometry for groups of input
  geometries. Use `""` to get a hull around all geometries combined. Use
  `NA` to get a hull around each geometry in `x`

- param:

  numeric. For `type="concave_ratio"`, use an edge length ratio between
  0 and 1. For `type="concave_length"`, the meaning depends on the input
  geometry: use a ratio between 0 and 1 for polygons, or a maximum edge
  length in coordinate units (0 or greater) for points and lines

- allowHoles:

  logical. May the output polygons contain holes? For "concave\_\*"
  methods only

- tight:

  logical. Should the hull follow the outer boundaries of the input
  polygons? For "concave_length" with polygon geometry only

- percent:

  numeric \> 0 and \<= 100. If below 100, only that percentage of
  geometries closest to the centroid is used

## Value

SpatVector

## Details

A concave hull is a polygon which contains all the points of the input.
It can be a better representation of the input data (typically points)
than the convex hull. There are many possible concave hulls with
different degrees of concaveness. These can be created with argument
`param`.

For `type="concave_ratio"`, larger ratios retain longer edges. A ratio
of 0 gives the most concave hull and 1 gives the convex hull. This
method requires GEOS 3.11 or later.

For `type="concave_length"` with polygons, the ratio applies to the
edges between polygons. At 0 the original polygons are retained. The
`tight` argument controls whether the hull follows their outer
boundaries. This method requires GEOS 3.11 or later. For points and
lines, a larger maximum edge length gives a less concave hull. This
method requires GEOS 3.12 or later.

## Examples

``` r
p <- vect(system.file("ex/lux.shp", package="terra"))
h <- hull(p)

plot(p)
lines(h, col="orange")

hh <- hull(p, "convex", by="NAME_1")
lines(hh, col="purple")


pts <- centroids(p)
plot(pts, ext=ext(p)+0.1)
lines(hull(pts, type="convex"), col="darkgreen")
lines(hull(pts, type="rectangle"), col="blue")
lines(hull(pts, type="circle"), col="red")
lines(hull(pts, percent=60), col="orange")

```
