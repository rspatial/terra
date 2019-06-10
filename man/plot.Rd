\name{plot}

\docType{methods}

\alias{plot}
\alias{plot,SpatRaster,missing-method}
\alias{plot,SpatVector,missing-method}
\alias{plot,SpatRaster,numeric-method}
\alias{image}
\alias{image,SpatRaster-method}
\alias{points,SpatVector-method}
\alias{lines,SpatVector-method}
\alias{lines,SpatExtent-method}
\alias{lines,SpatRaster-method}



\title{Plot a SpatRaster or SpatVector}

\description{
Plot (that is, make a map of) the values of a SpatRaster or SpatVector
}

\usage{
\S4method{plot}{SpatRaster,numeric}(x, y, cols, maxpixels=100000, leg.mar=3, 
		leg.levels=5, leg.shrink=c(0,0), leg.main=NULL, leg.main.cex=1,
		digits=1, useRaster=TRUE, zlim, xlab="", ylab="", ...)
\S4method{plot}{SpatRaster,missing}(x, y, maxpixels=100000, nc, nr, main, maxnl=16, ...)
\S4method{plot}{SpatVector,missing}(x, y, xlab="", ylab="", ...)
\S4method{image}{SpatRaster}(x, y, maxpixels=100000, xlab="", ylab="", useRaster=TRUE, ...)
\S4method{points}{SpatVector}(x, ...)
\S4method{lines}{SpatVector}(x, ...)
\S4method{lines}{SpatExtent}(x, ...)
}

\arguments{
  \item{x}{SpatRaster or SpatVector}
  \item{y}{missing or integer indicating the layer(s) to be plotted}
  \item{cols}{character. colors}
  \item{maxpixels}{integer > 0. Maximum number of cells to use for the plot}
  \item{leg.mar}{integer > 0. size of margin for the legend}
  \item{leg.levels}{integer > 0. number of labels on the legend}
  \item{leg.shrink}{numeric. relative shrinkage of legend from top and bottom}
  \item{leg.main}{character. legend title}
  \item{leg.main.cex}{character. legend title size}
  \item{digits}{integer >= 0}
  \item{zlim}{range of allowed values}
  \item{xlab}{character}
  \item{ylab}{character}
  \item{useRaster}{logical; if \code{TRUE} plotting is faster because a bitmap raster is used to plot the image instead of polygons. See \code{\link[graphics]{image}}}
  \item{nc}{integer > 0}
  \item{nr}{integer > 0}
  \item{main}{character}
  \item{maxnl}{integer > 0}
  \item{...}{additional graphical arguments}  
}


\examples{ 
f <- system.file("exdata/test.tif", package="terra") 
r <- rast(f)
plot(r)

f <- system.file("exdata/lux.shp", package="terra")
v <- vect(f)
r <- rast(v)
values(r) <- 1:ncell(r)
plot(r)
lines(v)
points(v)
}

\keyword{methods}
\keyword{spatial}

