\name{RGB}

\docType{methods}

\alias{RGB}
\alias{RGB,SpatRaster-method}
\alias{RGB<-}
\alias{RGB<-,SpatRaster-method}
\alias{colorize}
\alias{colorize,SpatRaster-method}
\alias{has.RGB}
\alias{has.RGB,SpatRaster-method}

  
\title{Layers representing colors}

\description{
With \code{RGB} you can get or set the layers to be used as Red, Green and Blue when plotting a SpatRaster. Currently, a benefit of this is that \code{\link[terra]{plot}} will send the object to \code{\link{plotRGB}}. You can also associated the layers with another color space (HSV, HSI or HSL)

With \code{colorize} you can convert a three-layer RGB SpatRaster into other color spaces. You can also convert it into a single-layer SpatRaster with a color-table.
}

\usage{
\S4method{RGB}{SpatRaster}(x)

\S4method{RGB}{SpatRaster}(x, ..., type="rgb")<-value

\S4method{colorize}{SpatRaster}(x, to="hsv", alpha=FALSE, stretch=NULL, 
	grays=FALSE, NAzero=FALSE, filename="", overwrite=FALSE, ...)

\S4method{has.RGB}{SpatRaster}(x, strict=TRUE)
}


\arguments{
  \item{x}{SpatRaster}
  
  \item{value}{vector of three (or four) positive integers indicating the layers that are red, green and blue (and optionally a fourth transparency layer). Or \code{NULL} to remove the RGB settings}
  
  \item{type}{character. The color space. One of "rgb" "hsv", "hsi" and "hsl"}

  \item{to}{character. The color space to transform the values to. If \code{x} has RGB set, you can transform these to "hsv", "hsi" and "hsl", or use "col" to create a single layer with a color table. You can also use "rgb" to back transform to RGB}

  \item{alpha}{logical. Should an alpha (transparency) channel be included? Only used if \code{x} has a color-table and \code{to="rgb"}}

  \item{stretch}{character. Option to stretch the values to increase contrast: "lin" (linear) or "hist" (histogram). Only used for transforming RGB to col}
  
  \item{grays}{logical. If \code{TRUE}, a gray-scale color-table is created. Only used for transforming RGB to col}

  \item{NAzero}{logical. If \code{TRUE}, NAs are treated as zeros such that a color can be returned if at least one of the three channels has a value. Only used for transforming RGB to (\code{"col"})}

  \item{strict}{logical. If \code{TRUE}, the function returns \code{FALSE} if a color space such as "hsv", "hsi" and "hsl" is used}
  
  \item{filename}{character. Output filename}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}  
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
}  


\examples{
r <- rast(system.file("ex/logo.tif", package="terra"))   
plot(r)
has.RGB(r)
RGB(r) <- NULL
has.RGB(r)
plot(r)
RGB(r) <- c(3,1,2)
plot(r)

RGB(r) <- 1:3
x <- colorize(r, "col")
y <- colorize(r, "hsv")
z <- colorize(y, "rgb")
}

\keyword{methods}
\keyword{spatial}

