\name{coerce}

\docType{methods}

\alias{as.vector}
\alias{as.matrix}
\alias{as.array}

\alias{as.vector,SpatRaster-method}
\alias{as.matrix,SpatRaster-method}
\alias{as.array,SpatRaster-method}


\title{Coercion of a SpatRaster to other object types}

\description{
Coercion to other object types
}

\usage{
\S4method{as.vector}{SpatRaster}(x, mode='any')

\S4method{as.matrix}{SpatRaster}(x, wide=FALSE)

\S4method{as.array}{SpatRaster}(x)
}


\arguments{
\item{x}{SpatRaster or SpatVector}
\item{mode}{this argument is ignored}
\item{wide}{logical}
}

\value{
vector, matrix, array
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

}

\keyword{spatial}

