\name{width}

\docType{methods}

\alias{width}
\alias{clearance}

\alias{width,SpatVector-method}
\alias{clearance,SpatVector-method}

\title{ 
SpatVector geometric properties
}

\description{
\code{width} returns the minimum diameter of the geometry, defined as the smallest band that contains the geometry, where a band is a strip of the plane defined by two parallel lines. This can be thought of as the smallest hole that the geometry can be moved through, with a single rotation.

\code{clearance} returns the minimum clearance of a geometry. The minimum clearance is the smallest amount by which a vertex could be moved to produce an invalid polygon, a non-simple linestring, or a multipoint with repeated points. If a geometry has a minimum clearance of 'mc', it can be said that:

No two distinct vertices in the geometry are separated by less than 'mc'
No vertex is closer than 'mc' to a line segment of which it is not an endpoint.
If the minimum clearance cannot be defined for a geometry (such as with a single point, or a multipoint whose points are identical, NA is returned.
}

\usage{
\S4method{width}{SpatVector}(x, as.lines=FALSE)
\S4method{clearance}{SpatVector}(x, as.lines=FALSE)
}

\arguments{
  \item{x}{SpatVector of lines or polygons}
  \item{as.lines}{logical. If \code{TRUE} lines are returned that define the width or clearance}
} 

\value{
numeric or SpatVector
}

\seealso{
\code{\link{minRect}}
}

\examples{
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
width(v)
clearance(v)
}

\keyword{methods}
\keyword{spatial}

