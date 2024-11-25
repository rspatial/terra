\name{rowSums}

\docType{methods}

\alias{rowSums}
\alias{rowSums,SpatRaster-method}
\alias{colSums}
\alias{colSums,SpatRaster-method}
\alias{rowMeans}
\alias{rowMeans,SpatRaster-method}
\alias{colMeans}
\alias{colMeans,SpatRaster-method}

\title{row/col sums and means for SpatRaster}

\description{
Sum or average values of SpatRaster layers by row or column. 
}

\usage{
\S4method{rowSums}{SpatRaster}(x, na.rm=FALSE, dims=1L, ...) 
\S4method{colSums}{SpatRaster}(x, na.rm=FALSE, dims=1L, ...) 
\S4method{rowMeans}{SpatRaster}(x, na.rm=FALSE, dims=1L, ...) 
\S4method{colMeans}{SpatRaster}(x, na.rm=FALSE, dims=1L, ...) 
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
r <- rast(ncols=2, nrows=5, nl=2, vals=1:20)
rowSums(r)
colSums(r)
colMeans(r)
}

\keyword{spatial}

