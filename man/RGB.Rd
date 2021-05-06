\name{RGB}

\docType{methods}

\alias{RGB}
\alias{RGB,SpatRaster-method}
\alias{RGB<-}
\alias{RGB<-,SpatRaster-method}

\alias{RGB2col}
\alias{RGB2col,SpatRaster-method}
  
\title{Declare RGB channels}

\description{
With \code{RGB} you can get or set the layers to be used as Red, Green and Blue when plotting a SpatRaster. Currently, the sole benefit of this is that \code{\link{plot}} will send the object to \code{\link{plotRGB}}

With \code{RGB2col} you can convert a three-layer RGB SpatRaster into a single-layer SpatRaster with a color-table.
}

\usage{
\S4method{RGB}{SpatRaster}(x)

\S4method{RGB}{SpatRaster}(x)<-value

\S4method{RGB2col}{SpatRaster}(x, value, stretch=NULL, grays=FALSE, filename="", overwrite=FALSE, ...)
}


\arguments{
  \item{x}{SpatRaster}
  \item{value}{vector of three positive integers indicating the layers that are red, green and blue. Or \code{NULL} to remove the RGB settings. With \code{rgb2coltab}, this can be missing if \code{RGB} is not \code{NULL
  \item{stretch}{character. Option to stretch the values to increase contrast: "lin" (linear) or "hist" (histogram)}  
  \item{grays}{logical. If \code{TRUE}, a gray-scale color-table is created}
  \item{filename}{character. Output filename}
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
}  


\examples{
s <- rast(system.file("ex/logo.tif", package="terra"))   
RGB(s) <- c(1,2,3)
plot(s)
RGB(s) <- NULL
}

\keyword{methods}
\keyword{spatial}

