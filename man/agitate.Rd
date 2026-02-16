\name{agitate}

\docType{methods}

\alias{agitate}
\alias{agitate,SpatVector-method}
\alias{agitate,ANY-method}

\title{Add noise to (jitter) a SpatVector of points}

\description{
Add noise to (jitter) a SpatVector of points. Points are moved randomly within a specified distance.
}

\usage{
\S4method{agitate}{SpatVector}(x, maxdist, ...)
}

\arguments{
  \item{x}{SpatVector. All other objects are passed to \code{\link{jitter}}}
  \item{maxdist}{positive number. The maximum distance from the original location, in meter for lon/lat, otherwise in the units of the CRS}
}

\value{
SpatVector
}

\examples{
## SpatVector
f <- system.file("ex/lux.shp", package="terra")
v <- as.points(vect(f)[1])
p <- agitate(v, 2500)
}

\keyword{spatial}
