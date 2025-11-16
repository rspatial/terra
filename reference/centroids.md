# Centroids

Get the centroids of polygons or lines, or centroid-like points that are
guaranteed to be inside the polygons or on the lines.

## Usage

``` r
# S4 method for class 'SpatVector'
centroids(x, inside=FALSE)
```

## Arguments

- x:

  SpatVector

- inside:

  logical. If `TRUE` the points returned are guaranteed to be inside the
  polygons or on the lines, but they are not the true centroids. True
  centroids may be outside a polygon, for example when a polygon is
  "bean shaped", and they are unlikely to be on their line

## Value

SpatVector of points

## Examples

``` r
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
x <- centroids(v)
y <- centroids(v, TRUE)
```
