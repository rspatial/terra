\name{rowSums}

\docType{methods}

\alias{rowSums}
\alias{rowSums,SpatRaster-method}
\alias{colSums}
\alias{colSums,SpatRaster-method}

\title{rowSums and colSums for SpatRaster}

\description{
Sum values of a SpatRaster by row or column. 
}

\usage{
\S4method{rowSums}{SpatRaster}(x, na.rm=FALSE, dims=1L,...) 
\S4method{colSums}{SpatRaster}(x, na.rm=FALSE, dims=1L,...) 
}

\arguments{
  \item{x}{SpatRaster}
  \item{na.rm}{logical. If \code{TRUE}, \code{NA} values are ignored}
  \item{dims}{this argument is ignored}
  \item{...}{additional arguments (none implemented)}
}

\value{
matrix
}

\seealso{
See \code{\link{global}} for summing all cells values
}

\examples{
r <- rast(ncols=2, nrows=5)
values(r) <- 1:10
as.matrix(r, wide=TRUE)
rowSums(r)
colSums(r)
}

\keyword{spatial}

