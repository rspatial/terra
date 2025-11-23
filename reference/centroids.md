# Centroids

Get the centroids of polygons or lines, or centroid-like points that are
guaranteed to be inside the polygons or on the lines.

Or get the (weighted) centroid of the the cells with values (not `NA`)
of a SpatRaster.

## Usage

``` r
# S4 method for class 'SpatVector'
centroids(x, inside=FALSE)

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
  "bean shaped", and they are unlikely to be on their line

- weighted:

  logical. If `TRUE` the centroids are computed as the weighted means of
  the coordinates of cells with values

SpatVector of points

f \<- system.file("ex/lux.shp", package="terra") v \<- vect(f) x \<-
centroids(v) y \<- centroids(v, TRUE)f \<- system.file("ex/elev.tif",
package="terra") r \<- rast(f) centroids(r)

spatial
