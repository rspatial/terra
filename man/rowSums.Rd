\name{rcSumMean}

\docType{methods}

\alias{rowSums}
\alias{rowMeans}
\alias{colSums}
\alias{colMeans}
\alias{rowSums,SpatRaster-method}
\alias{colSums,SpatRaster-method}
\alias{rowMeans,SpatRaster-method}
\alias{colMeans,SpatRaster-method}

\title{rowSums, colSums, rowMeans, colMeans for SpatRaster}

\description{
Sum or average values of a SpatRaster by row or column. 
}

\usage{
\S4method{rowSums}{SpatRaster}(x, na.rm=FALSE, dims=1L,...) 
\S4method{colSums}{SpatRaster}(x, na.rm=FALSE, dims=1L,...) 
\S4method{rowMeans}{SpatRaster}(x, na.rm=FALSE, dims=1L,...) 
\S4method{colMeans}{SpatRaster}(x, na.rm=FALSE, dims=1L,...) 
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
See \code{\link{global}} for summing or averaging all cell values
}

\examples{
r <- rast(ncols=2, nrows=5)
values(r) <- 1:10
as.matrix(r, wide=TRUE)
rowSums(r)
colSums(r)
rowMeans(r)
colMeans(r)
}

\keyword{spatial}

