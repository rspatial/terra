\docType{methods}
\name{watershed}

\alias{watershed}
\alias{watershed,SpatRaster-method}


\title{Catchment delineation}

\description{
delineate the area covered by a catchment from a SpatRaster with flow direction and a pour-point (catchment outlet).
}

\usage{
\S4method{watershed}{SpatRaster}(x, pourpoint, filename="",...) 
}

\arguments{
  \item{p}{SpatRaster with flow direction. See \code{\link{terrain}}. }
  \item{pourpoint}{matrix or SpatVector with the pour point location}
  \item{filename}{character. Output filename}
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
}

\value{ 
SpatRaster 
}

\author{
Ezio Crestaz, Emanuele Cordano, Roman Seliger
}

\examples{
elev <- rast(system.file('ex/elev_vinschgau.tif', package="terra"))
flowdir <- terrain(elev, "flowdir")
## pour point at Naturns 
pp <- cbind(653358.3, 5168222)
w <- watershed(flowdir, pp)
}

\keyword{spatial}
