\name{plotRGB}

\docType{methods}

\alias{plotRGB}
\alias{plotRGB,SpatRaster-method}

\title{Red-Green-Blue plot of a multi-layered SpatRaster}

\description{
Make a Red-Green-Blue plot based on three layers in a SpatRaster. The layers (sometimes referred to as "bands" because they may represent different bandwidths in the electromagnetic spectrum) are combined such that they represent the red, green and blue channel. This function can be used to make "true" (or "false") color images from Landsat and other multi-spectral satellite images. 

Note that the margins of the plot are set to zero (no axes or titles are visible) but can be set with the \code{mar} argument. 

An alternative way to plot RGB images is to first use \code{\link{colorize}} to create a single layer SpatRaster with a color-table and then use \code{\link[terra]{plot}}.

}
\usage{
\S4method{plotRGB}{SpatRaster}(x, r=1, g=2, b=3, a=NULL, scale=NULL, mar=0, 
		stretch=NULL, smooth=TRUE, colNA="white", alpha=NULL, bgalpha=NULL, 
		zlim=NULL, zcol=FALSE, axes=FALSE ,...)
}


\arguments{
  \item{x}{SpatRaster}
  \item{r}{integer between 1 and \code{nlyr(x)}. Layer to use as the Red channel}
  \item{g}{integer between 1 and \code{nlyr(x)}. Layer to use as the Green channel}
  \item{b}{integer between 1 and \code{nlyr(x)}. Layer to use as the Blue channel}
  \item{a}{NULL or integer between 1 and \code{nlyr(x)}. Layer to use as the alpha (transparency) channel. If not NULL, argument \code{alpha} is ignored}
  \item{scale}{integer. Maximum (possible) value in the three channels. Defaults to 255 or to the maximum value of \code{x} if that is known and larger than 255}
  \item{mar}{numeric vector recycled to length 4 to set the margins of the plot. Use \code{mar=NULL} or \code{mar=NA} to not set the margins}
  \item{stretch}{character. Option to stretch the values to increase contrast: "lin" (linear) or "hist" (histogram)}
  \item{smooth}{logical. If \code{TRUE}, smooth the image when drawing to get the appearance of a higher spatial resolution}
  \item{colNA}{color. The color used for cells that have NA values}
  \item{alpha}{transparency. Integer between 0 (transparent) and 255 (opaque)}
  \item{bgalpha}{Background transparency. Integer between 0 (transparent) and 255 (opaque)}
  \item{zlim}{numeric vector of length 2. Range of values to plot (optional). If this is set, and \code{stretch="lin"} is used, then the values are stretched within the range of \code{zlim}. This allows creating consistent coloring between SpatRasters with different cell-value ranges, even when stretching the colors for improved contrast}
  \item{zcol}{logical. If \code{TRUE} the values outside the range of zlim get the color of the extremes of the range. Otherwise, the values outside the zlim range get the color of \code{NA} values (see argument "colNA")}
  \item{axes}{logical. If \code{TRUE} axes are drawn (and arguments such as \code{main="title"} will be honored)}
  \item{...}{graphical parameters as in \code{\link{plot}}<SpatRaster-method>}
}


\seealso{ \code{\link{plot}}, \code{\link{colorize}}, \code{\link{RGB}} }

\examples{
b <- rast(system.file("ex/logo.tif", package="terra"))   
plotRGB(b)
plotRGB(b, mar=2)
plotRGB(b, 3, 2, 1)

b[1000:2000] <- NA
plotRGB(b, 3, 2, 1, stretch="hist")
}

\keyword{methods}
\keyword{spatial}

