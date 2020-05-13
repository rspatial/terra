\name{coerce}

\docType{methods}

\alias{as.vector}
\alias{as.matrix}
\alias{as.array}
\alias{as.polygons}
\alias{as.lines}
\alias{as.points}

\alias{as.vector,SpatRaster-method}
\alias{as.matrix,SpatRaster-method}
\alias{as.array,SpatRaster-method}
\alias{as.data.frame,SpatRaster-method}
\alias{as.polygons,SpatRaster-method}
\alias{as.points,SpatRaster-method}
\alias{as.lines,SpatVector-method}
\alias{as.points,SpatVector-method}


\title{Coercion to other object types}

\description{
Coercion to other object types or other vector types
}

\usage{
\S4method{as.vector}{SpatRaster}(x, mode='any')

\S4method{as.matrix}{SpatRaster}(x, wide=FALSE, ...)

\S4method{as.data.frame}{SpatRaster}(x, xy=FALSE, cells=FALSE, ...)

\S4method{as.array}{SpatRaster}(x, ...)

\S4method{as.polygons}{SpatRaster}(x, values=FALSE, na.rm=FALSE, ...)

\S4method{as.points}{SpatRaster}(x, values=FALSE, na.rm=FALSE, ...)

\S4method{as.lines}{SpatVector}(x, ...)

\S4method{as.points}{SpatVector}(x, ...)
}


\arguments{
\item{x}{SpatRaster or SpatVector}
\item{wide}{logical}
\item{xy}{logical}
\item{cells}{logical}
\item{mode}{this argument is ignored}
\item{values}{logical; include cell values as attributes?}
\item{na.rm}{logical; exclude cells with NAs?}
\item{...}{additional arguments. None implemented}
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
as.points(r)
p <- as.polygons(r)
p
as.lines(p)
as.points(p)
}

\keyword{spatial}

