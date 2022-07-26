\name{centroids}

\alias{centroids}
\alias{centroids,SpatVector-method}


\title{Centroids}

\description{
Get the centroids of polygons or lines, or centroid-like points that are guaranteed to be inside the polygons or on the lines.
}

\usage{
\S4method{centroids}{SpatVector}(x, inside=FALSE)
}

\arguments{
\item{x}{SpatVector}
\item{inside}{logical. If \code{TRUE} the points returned are guaranteed to be inside the polygons or on the lines, but they are not the true centroids. True centroids may be outside a polygon, for example when a polygon is "bean shaped", and they are unlikely to be on their line}
}


\value{SpatVector of points}


\examples{
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
x <- centroids(v)
y <- centroids(v, TRUE)
}

\keyword{spatial}
