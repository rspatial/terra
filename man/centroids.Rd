\name{centroids}

\alias{centroids}
\alias{centroids,SpatVector-method}
\alias{centroids,SpatRaster-method}


\title{Centroids}

\description{
SpatVector: Get the centroids of polygons or lines, or centroid-like points that are guaranteed to be inside the polygons or on the lines.

SpatRaster: Get the centroid of all cells, or he (weighted) centroid of the cells that are not \code{NA}.
}

\usage{
\S4method{centroids}{SpatVector}(x, inside=FALSE, correct=FALSE, geoms=TRUE, lonlat=FALSE, data.frame=FALSE)

\S4method{centroids}{SpatRaster}(x, weighted=FALSE, na.rm=TRUE, lonlat=FALSE, data.frame=FALSE)
}

\arguments{
\item{x}{SpatVector}
\item{inside}{logical. If \code{TRUE} the points returned are guaranteed to be inside the polygons or on the lines, but they are not the true centroids. True centroids may be outside a polygon, for example when a polygon is "bean shaped".}
\item{correct}{If \code{TRUE} corrected centroids are returned for geometries in \code{x} that have a true centroid that is not inside the geometry. If \code{inside=TRUE} it is replaced with another point that is guaranteed to be inside the geometry. If it is \code{FALSE}, it is replaced with the nearest point on the border of the geometry}
\item{geoms}{If \code{TRUE} centroids are returned for each geometry in \code{x}. Otherwise, a single centroid is returned for the aggregated geometries}
\item{weighted}{logical. If \code{TRUE} the centroids are computed as the weighted means of the coordinates of cells with values}
\item{na.rm}{logical. If \code{TRUE} cells that are \code{NA} are not included in the computation. Cannot be combined with \code{na.rm=FALSE}}
\item{lonlat}{logical. If \code{TRUE} the centroids are always longitude/latitude coordinates, irrespective of the coordinate reference system of \code{x}}
\item{data.frame}{logical. If \code{TRUE} a \code{data.frame} is returned instead of a \code{SpatVector}}
}

\value{SpatVector of points or data.frame}


\examples{
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
x <- centroids(v)
y <- centroids(v, TRUE)
z <- centroids(v, correct=TRUE) ## same as 'x' in this case


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
