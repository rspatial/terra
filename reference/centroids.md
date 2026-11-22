# Centroids

SpatVector: Get the centroids of polygons or lines, or centroid-like
points that are guaranteed to be inside the polygons or on the lines.

SpatRaster: Get the centroid of all cells, or he (weighted) centroid of
the cells that are not `NA`.

## Usage

``` r
# S4 method for class 'SpatVector'
centroids(x, inside=FALSE, correct=FALSE, geoms=TRUE, lonlat=FALSE, data.frame=FALSE)

# S4 method for class 'SpatRaster'
centroids(x, weighted=FALSE, na.rm=TRUE, lonlat=FALSE, data.frame=FALSE)
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

- geoms:

  If `TRUE` centroids are returned for each geometry in `x`. Otherwise,
  a single centroid is returned for the aggregated geometries

- weighted:

  logical. If `TRUE` the centroids are computed as the weighted means of
  the coordinates of cells with values

- na.rm:

  logical. If `TRUE` cells that are `NA` are not included in the
  computation. Cannot be combined with `na.rm=FALSE`

- lonlat:

  logical. If `TRUE` the centroids are always longitude/latitude
  coordinates, irrespective of the coordinate reference system of `x`

- data.frame:

  logical. If `TRUE` a `data.frame` is returned instead of a
  `SpatVector`

## Value

SpatVector of points or data.frame

## Examples

``` r
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
x <- centroids(v)
y <- centroids(v, TRUE)
z <- centroids(v, correct=TRUE) ## same as 'x' in this case


f <- system.file("ex/elev.tif", package="terra")
r <- rast(f)
centroids(r)
#> class       : SpatVector
#> geometry    : points
#> dimensions  : 1, 0  (geometries, attributes)
#> extent      : 6.091757, 6.091757, 49.77608, 49.77608  (xmin, xmax, ymin, ymax)
#> coord. ref. : 


p <- vect("POLYGON ((5.9 50, 6.3 50, 6.5 49.8, 6.3 49.6, 5.9 49.5, 6.3 49.7, 6.3 49.9, 5.9 50))")
plot(p)
points(centroids(p), col="red", cex=2)
points(centroids(p, inside=TRUE), col="blue", cex=2)
points(centroids(p, inside=FALSE, correct=TRUE), col="green", cex=2)
```
