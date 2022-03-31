\name{centroids}

\alias{centroids}
\alias{centroids,SpatVector-method}


\title{Centroids}

\description{
Get the centroids of polygons, or centroid-like points that are guaranteed to be inside the polygons
}

\usage{
\S4method{centroids}{SpatVector}(x, inside=FALSE)
}

\arguments{
\item{x}{SpatVector}
\item{inside}{logical. If \code{TRUE} the points returned are not the true centroids, but they are guaranteed to be inside the polygons. True centroids that may be outside a polygon (e.g. when a polygon is "bean shaped"}
}


\value{SpatVector of points}


\examples{
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
x <- centroids(v)
y <- centroids(v, TRUE)
}

\keyword{spatial}
