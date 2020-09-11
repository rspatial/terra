\name{c}

\docType{methods}


\alias{c}
\alias{c,SpatRaster-method}
\alias{c,SpatDataSet-method}

\alias{add<-}
\alias{add<-,SpatRaster-method}

\alias{collapse,SpatRaster-method}

\title{Combine SpatRasters}

\description{
Combine SpatRaster objects into a new object with \code{c} or add to an existing object with \code{add}. \code{collapse} combines sources within an object (that are in memory, or from the same file) to allow for faster processing.
}

\usage{
\S4method{c}{SpatRaster}(x, ...)
\S4method{add}{SpatRaster}(x)<-value
\S4method{collapse}{SpatRaster}(x, ...)

\S4method{c}{SpatDataSet}(x, ...)
}

\arguments{
  \item{x}{SpatRaster or SpatDataSet}
  \item{...}{For \code{c}: SpatRaster objects to be combined with \code{x}. Not implemented for \code{collapse)}
}



\value{
SpatRaster
}


\examples{
r <- rast(nrow=5, ncol=9)
values(r) <- 1:ncell(r)
x <- c(r, r*2, r*3)
}

\keyword{spatial}
