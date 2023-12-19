\name{all.equal}

\docType{methods}

\alias{all.equal}
\alias{all.equal,SpatRaster,SpatRaster-method}

\title{Compare two SpatRasters for equality}

\description{
Compare two SpatRasters for (near) equality. 

First the attributes of the objects are compared. If these are the same, a (perhaps small) sample of the raster cells is compared as well. 

The sample size used can be increased with the \code{maxcell} argument. You can set it to \code{Inf}, but for large rasters your computer may not have sufficient memory. See the examples for a safe way to compare all values. 
}

\usage{
\S4method{all.equal}{SpatRaster,SpatRaster}(target, current, maxcell=100000, ...)
}

\arguments{
  \item{target}{SpatRaster}
  \item{current}{SpatRaster}
  \item{maxcell}{positive integer. The size of the regular sample used to compare cell values}
  \item{...}{additional arguments passed to \code{\link{all.equal.numeric}} to compare cell values}
}

\seealso{\code{\link{identical}}, \code{\link{compareGeom}}}

\value{
Either \code{TRUE} or a character vector describing the differences between target and current.
}

\examples{
x <- sqrt(1:100)
mat <- matrix(x, 10, 10)
r1 <- rast(nrows=10, ncols=10, xmin=0, vals = x)
r2 <- rast(nrows=10, ncols=10, xmin=0, vals = mat)

all.equal(r1, r2)
all.equal(r1, r1*1)
all.equal(rast(r1), rast(r2))

# compare geometries 
compareGeom(r1, r2)

# Compare all cell values for near equality
# as floating point number imprecision can be a problem
m <- minmax(r1 - r2)
all(abs(m) < 1e-7)

# comparison of cell values to create new SpatRaster
e <- r1 == r2
}

\keyword{spatial}
