\name{intersect}

\docType{methods}

\alias{intersect}
\alias{intersects}

\alias{intersect,SpatVector,SpatVector-method}
\alias{intersect,SpatVector,SpatExtent-method}
\alias{intersect,SpatExtent,SpatVector-method}
\alias{intersect,SpatExtent,SpatExtent-method}
\alias{intersects,SpatVector,SpatVector-method}

\title{ 
Intersection
}

\description{
Intersect SpatVector, typically polygons or, with \code{intersects}, get a matrix indicating which geometries intersect.

Intersecting points with points uses the extent of \code{y} to get the intersection. Intersecting of points and lines is not supported because of numerical inaccuracies with that. You can use \code{\link{buffer}}, to create polygons from lines and use these with intersect. 

See \code{\link{crop}} for intersection of a SpatRaster.
}

\usage{
\S4method{intersect}{SpatVector,SpatVector}(x, y)

\S4method{intersect}{SpatVector,SpatExtent}(x, y)

\S4method{intersect}{SpatExtent,SpatVector}(x, y)

\S4method{intersect}{SpatExtent,SpatExtent}(x, y)

\S4method{intersects}{SpatVector,SpatVector}(x, y, ...)
}

\arguments{
  \item{x}{SpatVector or SpatExtent}
  \item{y}{SpatVector or SpatExtent}
} 

\value{
Same as \code{x}
}

\seealso{
\code{\link{union}}, code{\link{crop}}
}

\examples{
e1 <- ext(-10, 10, -20, 20)
e2 <- ext(0, 20, -40, 5)
intersect(e1, e2)
}


\keyword{methods}
\keyword{spatial}

