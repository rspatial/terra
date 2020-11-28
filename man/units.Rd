\name{units}

\alias{units}
\alias{units<-}
\alias{units,SpatRaster-method}
\alias{units<-,SpatRaster-method}

  
\title{units of SpatRaster layers}

\description{
Get or set the units of the layers of a SpatRaster.
}

\usage{
\S4method{units}{SpatRaster}(x, ...)

\S4method{units}{SpatRaster}(x)<-value

}

\arguments{
  \item{x}{SpatRaster}
  \item{value}{character vector}
  \item{...}{additional arguments. None implemented}
}

\value{
character  
}

\seealso{\code{\link{time}, \link{names}}}

\examples{
s <- rast(system.file("ex/logo.tif", package="terra"))   

units(s) <- rep("m/s", nlyr(s))
units(s)
s
}

\keyword{spatial}

