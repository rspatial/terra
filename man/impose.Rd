\name{impose}

\docType{methods}

\alias{impose}
\alias{impose,SpatRasterCollection-method}

\title{ 
Impose the geometry of a SpatRaster to those in a SpatRasterCollection.
}

\description{
Add empty layers in bewteen existing layers such that the time step between each layer is the same. 
See \code{\link{approximate}} to estimate values for these layer (and other missing values)
}

\usage{
\S4method{impose}{SpatRasterCollection}(x, y, filename="", ...)
}

\arguments{
  \item{x}{SpatRasterCollection}
  \item{y}{SpatRaster}
  \item{filename}{character. Output filename}
  \item{...}{list with named options for writing files as in \code{\link{writeRaster}}}
} 

\value{
SpatRaster
}

\seealso{
\code{\link{resample}}
}


\keyword{methods}
\keyword{spatial}

