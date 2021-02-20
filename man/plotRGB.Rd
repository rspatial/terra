\name{plotRGB}

\docType{methods}

\alias{plotRGB}
\alias{plotRGB,SpatRaster-method}

\title{Red-Green-Blue plot of a multi-layered SpatRaster}

\description{
Make a Red-Green-Blue plot based on three layers in a SpatRaster. The layers (sometimes referred to as "bands" because they may represent different bandwidths in the electromagnetic spectrum) are combined such that they represent the red, green and blue channel. This function can be used to make "true" (or "false") color images from Landsat and other multi-spectral satellite images. 
}



\usage{
\S4method{plotRGB}{SpatRaster}(x, r=1, g=2, b=3, scale, maxcell=500000, mar=0, stretch=NULL, 
ext=NULL, interpolate=FALSE, colNA="white", alpha, bgalpha, addfun=NULL, zlim=NULL, 
zlimcol=NULL, axes=FALSE, xlab="", ylab="", asp=NULL, add=FALSE, ...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{r}{integer. Index of the Red channel, between 1 and nlyr(x)}
  \item{g}{integer. Index of the Green channel, between 1 and nlyr(x)}
  \item{b}{integer. Index of the Blue channel, between 1 and nlyr(x)}
  \item{scale}{integer. Maximum (possible) value in the three channels. Defaults to 255 or to the maximum value of \code{x} if that is known and larger than 255}
  \item{maxcell}{integer > 0. Maximum number of pixels to use}
  \item{mar}{numeric vector recycled to length 4 to set the margins of the plot. Use \code{mar=NULL} or \code{mar=NA} to not set the margins}
  \item{stretch}{character. Option to stretch the values to increase the contrast of the image: "lin" or "hist"}
  \item{ext}{An \code{\link{SpatExtent}} object to zoom in to a region of interest (see \code{\link{draw}})}
  \item{interpolate}{logical. If \code{TRUE}, interpolate the image when drawing}
  \item{colNA}{color for the background (\code{NA} values)}
  \item{alpha}{transparency. Integer between 0 (transparent) and 255 (opaque)}
  \item{bgalpha}{Background transparency. Integer between 0 (transparent) and 255 (opaque)}
  \item{addfun}{Function to add additional items such as points or polygons to the plot (map). See \code{\link{plot}}}
  \item{zlim}{numeric vector of length 2. Range of values to plot (optional)}
  \item{zlimcol}{If \code{NULL} the values outside the range of zlim get the color of the extremes of the range. If zlimcol has any other value, the values outside the zlim range get the color of \code{NA} values (see colNA)}
  \item{axes}{logical. If \code{TRUE} axes are drawn (and arguments such as \code{main="title"} will be honored)}
  \item{xlab}{character. Label of x-axis}
  \item{ylab}{character. Label of y-axis}
  \item{asp}{numeric. Aspect (ratio of x and y. If NULL, and appropriate value is computed to match data for the longitude/latitude coordinate reference system, and 1 for planar coordinate reference systems}
  \item{add}{logical. If \code{TRUE} add values to current plot}
  \item{margins}{logical. If \code{TRUE} standard whitespace margins are used. If \code{FALSE}, graphics::par(plt=c(0,1,0,1)) is used}
  \item{...}{graphical parameters as in \code{\link{plot}} or \code{\link{rasterImage}}}
}


\seealso{ \code{\link{plot}} }

\examples{
b <- rast(system.file("ex/logo.tif", package="terra"))   
plotRGB(b)
graphics::par(plt=c(0,1,0,1))
plotRGB(b, 3, 2, 1)
plotRGB(b, 3, 2, 1, stretch='hist')
}

\keyword{methods}
\keyword{spatial}

