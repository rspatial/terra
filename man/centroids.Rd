\name{centroids}

\alias{centroids}
\alias{centroids,SpatVector-method}
\alias{centroids,SpatRaster-method}


\title{Centroids}

\description{
Get the centroids of polygons or lines, or centroid-like points that are guaranteed to be inside the polygons or on the lines.

Or get the (weighted) centroid of the the cells with values (not \code{NA}) of a SpatRaster.
}

\usage{
\S4method{centroids}{SpatVector}(x, inside=FALSE)

\S4method{centroids}{SpatRaster}(x, weighted=FALSE)
}

\arguments{
\item{x}{SpatVector}
\item{inside}{logical. If \code{TRUE} the points returned are guaranteed to be inside the polygons or on the lines, but they are not the true centroids. True centroids may be outside a polygon, for example when a polygon is "bean shaped", and they are unlikely to be on their line}
\item{weighted}{logical. If \code{TRUE} the centroids are computed as the weighted means of the coordinates of cells with values}
}

\value{SpatVector of points}


\examples{
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
x <- centroids(v)
y <- centroids(v, TRUE)


f <- system.file("ex/elev.tif", package="terra")
r <- rast(f)
centroids(r)
}

\keyword{spatial}
