\name{all.equal}

\docType{methods}

\alias{all.equal}
\alias{all.equal,SpatRaster,SpatRaster-method}

\title{Compare two SpatRasters for equality}

\description{
Compare two SpatRasters for (near) equality
}

\usage{
\S4method{all.equal}{SpatRaster,SpatRaster}(current, target, maxcell=10000)
}

\arguments{
  \item{target}{SpatRaster}
  \item{current}{SpatRaster}
  \item{maxcell}{postive integer. The sample size to compare cell values}
}

\seealso{\code{\link{compareGeom}}}

\value{
Either \code{TRUE} or a vector of mode "character" describing the differences between target and current.
}

\examples{
x <- 1:100
mat <- matrix(x, 10, 10)
r1 <- rast(nrows=10, ncols=10, xmin=0, vals = x)
r2 <- rast(nrows=10, ncols=10, xmin=0, vals = mat)

all.equal(r1, r2)
all.equal(r1, r1*1)
all.equal(rast(r1), rast(r2))
}

\keyword{spatial}
