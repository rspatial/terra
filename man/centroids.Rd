\name{centroids}

\alias{centroids}
\alias{centroids,SpatVector-method}
\alias{centroids,SpatRaster-method}


\title{Centroids}

\description{
Get the centroids of polygons or lines, or centroid-like points that are guaranteed to be inside the polygons or on the lines.

Or get the (weighted) centroid of the cells with values (not \code{NA}) of a SpatRaster.
}

\usage{
\S4method{centroids}{SpatVector}(x, inside=FALSE, correct=FALSE)

\S4method{centroids}{SpatRaster}(x, weighted=FALSE)
}

\arguments{
\item{x}{SpatVector}
\item{inside}{logical. If \code{TRUE} the points returned are guaranteed to be inside the polygons or on the lines, but they are not the true centroids. True centroids may be outside a polygon, for example when a polygon is "bean shaped".}
\item{correct}{If \code{TRUE} corrected centroids are returned for geometries in \code{x} that have a true centroid that is not inside the geometry. If \code{inside=TRUE} it is replaced with another point that is guaranteed to be inside the geometry. If it is \code{FALSE}, it is replaced with the nearest point on the border of the geometry}
\item{weighted}{logical. If \code{TRUE} the centroids are computed as the weighted means of the coordinates of cells with values}
}

\value{SpatVector of points}


\examples{
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
x <- centroids(v)
y <- centroids(v, TRUE)
z <- centroids(v, correct=TRUE) ## same as z in this case


f <- system.file("ex/elev.tif", package="terra")
r <- rast(f)
centroids(r)


p <- vect("POLYGON ((5.9 50, 6.3 50, 6.5 49.8, 6.3 49.6, 5.9 49.5, 6.3 49.7, 6.3 49.9, 5.9 50))")
plot(p)
points(centroids(p), col="red", cex=2)
points(centroids(p, inside=TRUE), col="blue", cex=2)
points(centroids(p, inside=FALSE, correct=TRUE), col="green", cex=2)
}

\keyword{spatial}