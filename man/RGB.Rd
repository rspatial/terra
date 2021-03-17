\name{RGB}

\docType{methods}

\alias{RGB<-}
\alias{RGB<-,SpatRaster-method}
  
\title{Declare RGB channels}

\description{
Declare the layers to be used as Red, Green and Blue when plotting a SpatRaster. Currently, the sole benefit of this is that \code{\link{plot}} will send the object to \code{\link{plotRGB}}
}

\usage{

\S4method{RGB}{SpatRaster}(x)<-value
}


\arguments{
  \item{x}{SpatRaster}
  \item{value}{positive integers indicating the layers that are red, green and blue. Or \code{NULL} to remove the RGB settings}
}  


\examples{
s <- rast(system.file("ex/logo.tif", package="terra"))   
RGB(s) <- c(1,2,3)
plot(s)
RGB(s) <- NULL
}

\keyword{methods}
\keyword{spatial}

