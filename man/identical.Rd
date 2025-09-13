\name{identical}

\docType{methods}

\alias{identical}
\alias{identical,SpatRaster,SpatRaster-method}
\alias{identical,SpatVector,SpatVector-method}
\alias{identical,SpatExtent,SpatExtent-method}

\title{Compare two SpatRaster, SpatVector or SpatExtent objects for equality}

\description{
When, comparing two SpatRasters for equality, first the attributes of the objects are compared. If these are the same, a the raster cells are compared as well. This can be time consuming, and you may prefer to use a sample instead with \code{\link{all.equal}}
}

\usage{
\S4method{identical}{SpatRaster,SpatRaster}(x, y)

\S4method{identical}{SpatVector,SpatVector}(x, y)

\S4method{identical}{SpatExtent,SpatExtent}(x, y)
}

\arguments{
  \item{x}{SpatRaster, SpatVector, or SpatExtent}
  \item{y}{object of the same class as \code{x}}
}

\seealso{\code{\link{all.equal}}, \code{\link{compareGeom}}}

\value{
single logical value
}

\examples{
x <- sqrt(1:100)
mat <- matrix(x, 10, 10)
r1 <- rast(nrows=10, ncols=10, xmin=0, vals = x)
r2 <- rast(nrows=10, ncols=10, xmin=0, vals = t(mat))

identical(r1, r2)
identical(r1, r1*1)
identical(rast(r1), rast(r2))
}

\keyword{spatial}
