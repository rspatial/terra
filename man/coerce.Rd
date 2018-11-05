\name{coerce}

\docType{methods}

\alias{as.vector}
\alias{as.matrix}
\alias{as.array}

\alias{as.vector,SpatRaster-method}
\alias{as.matrix,SpatRaster-method}
\alias{as.array,SpatRaster-method}
\alias{as.data.frame,SpatRaster-method}



\title{Coercion to other object types}

\description{
Coercion to other object types
}

\usage{
\S4method{as.vector}{SpatRaster}(x, mode='any')
\S4method{as.matrix}{SpatRaster}(x, wide=FALSE, ...)
\S4method{as.data.frame}{SpatRaster}(x, xy=FALSE, cells=FALSE, ...)
\S4method{as.array}{SpatRaster}(x, ...)
}


\arguments{
\item{x}{SpatRaster object}
\item{wide}{Logical}
\item{xy}{Logical}
\item{cells}{Logical}
\item{mode}{this argument is ignored}
\item{...}{Additional arguments. None implemented}
}

\value{
vector, matrix, array or data.frame
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

