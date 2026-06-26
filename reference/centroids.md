# Centroids

Get the centroids of polygons or lines, or centroid-like points that are
guaranteed to be inside the polygons or on the lines.

Or get the (weighted) centroid of the cells with values (not `NA`) of a
SpatRaster.

## Usage

``` r
# S4 method for class 'SpatVector'
centroids(x, inside=FALSE, correct=FALSE)

# S4 method for class 'SpatRaster'
centroids(x, weighted=FALSE)
```

## Arguments

- x:

  SpatVector

- inside:

  logical. If `TRUE` the points returned are guaranteed to be inside the
  polygons or on the lines, but they are not the true centroids. True
  centroids may be outside a polygon, for example when a polygon is
  "bean shaped".

- correct:

  If `TRUE` corrected centroids are returned for geometries in `x` that
  have a true centroid that is not inside the geometry. If `inside=TRUE`
  it is replaced with another point that is guaranteed to be inside the
  geometry. If it is `FALSE`, it is replaced with the nearest point on
  the border of the geometry

- weighted:

  logical. If `TRUE` the centroids are computed as the weighted means of
  the coordinates of cells with values

## Value

SpatVector of points

## Examples

``` r
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
x <- centroids(v)
y <- centroids(v, TRUE)
z <- centroids(v, correct=TRUE) ## same as z in this case


f <- system.file("ex/elev.tif", package="terra")
r <- rast(f)
centroids(r)
#>          x        y
#> 1 6.091757 49.77608


p <- vect("POLYGON ((5.9 50, 6.3 50, 6.5 49.8, 6.3 49.6, 5.9 49.5, 6.3 49.7, 6.3 49.9, 5.9 50))")
plot(p)
points(centroids(p), col="red", cex=2)
points(centroids(p, inside=TRUE), col="blue", cex=2)
points(centroids(p, inside=FALSE, correct=TRUE), col="green", cex=2)
```
