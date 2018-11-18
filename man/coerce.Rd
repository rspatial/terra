\name{coerce}

\docType{methods}

\alias{as.vector}
\alias{as.matrix}
\alias{as.array}
\alias{as.polygons}

\alias{as.vector,SpatRaster-method}
\alias{as.matrix,SpatRaster-method}
\alias{as.array,SpatRaster-method}
\alias{as.data.frame,SpatRaster-method}
\alias{as.polygons,SpatRaster-method}


\title{Coercion to other object types}

\description{
Coercion to other object types
}

\usage{
\S4method{as.vector}{SpatRaster}(x, mode='any')
\S4method{as.matrix}{SpatRaster}(x, wide=FALSE, ...)
\S4method{as.data.frame}{SpatRaster}(x, xy=FALSE, cells=FALSE, ...)
\S4method{as.array}{SpatRaster}(x, ...)
\S4method{as.polygons}{SpatRaster}(x, values=FALSE, na.rm=FALSE, ...)
}


\arguments{
\item{x}{SpatRaster}
\item{wide}{logical}
\item{xy}{logical}
\item{cells}{logical}
\item{mode}{this argument is ignored}
\item{values}{logical; include cell values as attributes?}
\item{na.rm}{logical; exclude cells with NAs?}
\item{...}{Additional arguments. None implemented}
}

\value{
vector, matrix, array, data.frame or SpatVector
}

\examples{
r <- rast(ncol=2, nrow=2)
values(r) <- 1:ncell(r)

as.vector(r)
as.matrix(r)
as.matrix(r, wide=TRUE)
as.data.frame(r, xy=TRUE)
as.array(r)

}

\keyword{spatial}

