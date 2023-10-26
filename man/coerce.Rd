\name{coerce}

\docType{methods}

\alias{as.vector}
\alias{as.matrix}
\alias{as.array}

\alias{as.vector,SpatRaster-method}
\alias{as.matrix,SpatRaster-method}
\alias{as.array,SpatRaster-method}

\alias{as.vector,SpatExtent-method}
\alias{as.matrix,SpatExtent-method}


\title{Coercion to vector, matrix or array}

\description{
Coercion of a SpatRaster to a vector, matrix or array. Or coerce a SpatExtent to a vector or matrix
}

\usage{
\S4method{as.vector}{SpatRaster}(x, mode='any')

\S4method{as.matrix}{SpatRaster}(x, wide=FALSE, ...)

\S4method{as.array}{SpatRaster}(x)

\S4method{as.vector}{SpatExtent}(x, mode='any')

\S4method{as.matrix}{SpatExtent}(x, ...)
}


\arguments{
\item{x}{SpatRaster or SpatVector}
\item{wide}{logical. If \code{FALSE} each layer in the SpatRaster becomes a column in the matrix and each cell in the SpatRaster becomes a row. If \code{TRUE} each row in the SpatRaster becomes a row in the matrix and each column in the SpatRaster becomes a column in the matrix }
\item{mode}{this argument is ignored}
\item{...}{additional arguments (none implemented)}
}

\value{
vector, matrix, or array
}

\seealso{\code{\link{as.data.frame}} and \code{\link{as.polygons}}}


\examples{
r <- rast(ncols=2, nrows=2)
values(r) <- 1:ncell(r)

as.vector(r)
as.matrix(r)
as.matrix(r, wide=TRUE)
as.data.frame(r, xy=TRUE)
as.array(r)

as.vector(ext(r))
as.matrix(ext(r))
}

\keyword{spatial}

