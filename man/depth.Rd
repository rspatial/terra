\name{depth}

\alias{depth}
\alias{depth<-}
\alias{depth,SpatRaster-method}
\alias{depth<-,SpatRaster-method}
\alias{depthName}
\alias{depthName<-}
\alias{depthName,SpatRaster-method}
\alias{depthName<-,SpatRaster-method}
\alias{depthUnit}
\alias{depthUnit<-}
\alias{depthUnit,SpatRaster-method}
\alias{depthUnit<-,SpatRaster-method}

  
\title{depth of SpatRaster layers}

\description{
Get or set the depth of the layers of a SpatRaster. Experimental. 
}

\usage{
\S4method{depth}{SpatRaster}(x)

\S4method{depth}{SpatRaster}(x)<-value

\S4method{depthName}{SpatRaster}(x)

\S4method{depthName}{SpatRaster}(x)<-value

\S4method{depthUnit}{SpatRaster}(x)

\S4method{depthUnit}{SpatRaster}(x)<-value
}

\arguments{
  \item{x}{SpatRaster}
  \item{value}{numeric vector (\code{depth}), or character (\code{depthName} and \code{depthUnit})}
}

\value{
numeric or SpatRaster
}

\seealso{\code{\link{time}}}

\examples{
s <- rast(nlyr=3)

depth(s) <- c(0, pi/2, pi)
depth(s)

depthName(s) <- "angle"
depthUnit(s) <- "radians"
s
}

\keyword{spatial}

