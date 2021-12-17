\name{RGB}

\docType{methods}

\alias{RGB}
\alias{RGB,SpatRaster-method}
\alias{RGB<-}
\alias{RGB<-,SpatRaster-method}
\alias{colorize}
\alias{colorize,SpatRaster-method}

\alias{RGB2col}
\alias{RGB2col,SpatRaster-method}
  
\title{Layers representing colors}

\description{
With \code{RGB} you can get or set the layers to be used as Red, Green and Blue when plotting a SpatRaster. Currently, a benefit of this is that \code{\link[terra]{plot}} will send the object to \code{\link{plotRGB}}

With \code{colorize} you can convert a three-layer RGB SpatRaster into other color spaces. You can also convert it into a single-layer SpatRaster with a color-table.
}

\usage{
\S4method{RGB}{SpatRaster}(x)

\S4method{RGB}{SpatRaster}(x)<-value

\S4method{colorize}{SpatRaster}(x, to, stretch=NULL, grays=FALSE, filename="", overwrite=FALSE, ...)
}


\arguments{
  \item{x}{SpatRaster}
  
  \item{value}{vector of three (or four) positive integers indicating the layers that are red, green and blue (and optionally a fourth transparancy layer). Or \code{NULL} to remove the RGB settings}
  
  \item{to}{character. The color space to transform the values to. If \code{x} has RGB set, you can transform these to "hsv", "hsi" and "hsl", or use "col" to create a single layer with a color table. You can also use "rgb" to backtransform to RGB}

  \item{stretch}{character. Option to stretch the values to increase contrast: "lin" (linear) or "hist" (histogram). Only used for transforming RGB to col}
  
  \item{grays}{logical. If \code{TRUE}, a gray-scale color-table is created. Only used for transforming RGB to col}
  
  \item{filename}{character. Output filename}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}  
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
}  


\examples{
s <- rast(system.file("ex/logo.tif", package="terra"))   
plot(s)
RGB(s) <- NULL
plot(s)
RGB(s) <- c(3,1,2)
plot(s)

RGB(s) <- 1:3
x <- colorize(x, "col")
y <- colorize(x, "hsv")
z <- colorize(y, "rgb")
}

\keyword{methods}
\keyword{spatial}

