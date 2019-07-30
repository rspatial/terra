\name{plot}

\docType{methods}

\alias{plot}
\alias{plot,SpatRaster,missing-method}
\alias{plot,SpatRaster,numeric-method}
\alias{plot,SpatRaster,SpatRaster-method}
\alias{plot,SpatVector,missing-method}
\alias{plot,SpatExtent,missing-method}

\alias{image}
\alias{image,SpatRaster-method}

\alias{lines,SpatRaster-method}
\alias{lines,SpatVector-method}
\alias{lines,SpatExtent-method}

\alias{points,SpatVector-method}
\alias{points,SpatExtent-method}


\title{Plot a SpatRaster or SpatVector}

\description{
Plot (that is, make a map of) the values of a SpatRaster or SpatVector
}

\usage{
\S4method{plot}{SpatRaster,numeric}(x, y, col, maxcell=100000, leg.mar, 
		leg.levels=5, leg.shrink=c(0,0), leg.main=NULL, leg.main.cex=1, leg.ext,
		digits, useRaster=TRUE, zlim, xlab="", ylab="", ...)

\S4method{plot}{SpatRaster,missing}(x, y, maxcell=100000, nc, nr, main, maxnl=16, ...)

\S4method{plot}{SpatRaster,SpatRaster}(x, y, maxcell=100000, nc, nr, maxnl=16, gridded=FALSE, ncol=25, nrow=25, ...)

\S4method{image}{SpatRaster}(x, y, maxcell=100000, xlab="", ylab="", ...)

\S4method{plot}{SpatVector,missing}(x, y, col=NULL, xlab="", ylab="", add=FALSE, ...)
\S4method{points}{SpatVector}(x, col=NULL, ...)
\S4method{lines}{SpatVector}(x, col=NULL, ...)

\S4method{plot}{SpatExtent,missing}(x, col=NULL, ...)
\S4method{lines}{SpatExtent}(x, col=NULL, ...)
\S4method{points}{SpatExtent}(x, col=NULL, ...)
}

\arguments{
  \item{x}{SpatRaster or SpatVector}
  \item{y}{missing or integer indicating the layer(s) to be plotted}
  \item{col}{character. colors}
  \item{maxcell}{integer > 0. Maximum number of cells to use for the plot}
  \item{leg.mar}{integer > 0. size of margin for the legend}
  \item{leg.levels}{integer > 0. number of labels on the legend}
  \item{leg.shrink}{numeric. relative shrinkage of legend from top and bottom}
  \item{leg.main}{character. legend title}
  \item{leg.main.cex}{character. legend title size}
  \item{leg.ext}{SpatExtent or numeric vector indicating the area where to plot the legend. The vector should have a lenght of 4 (xmin, xmax, ymin, ymax)}
  \item{digits}{integer >= 0}
  \item{zlim}{range of allowed values}
  \item{xlab}{character}
  \item{ylab}{character}
  \item{useRaster}{logical; if \code{TRUE} plotting is faster because a bitmap raster is used to plot the image instead of polygons. See \code{\link[graphics]{image}}}
  \item{nc}{integer > 0. Optional. The number of columns to divide the plotting device in (when plotting multiple layers)}
  \item{nr}{integer > 0. Optional. The number of rows to divide the plotting device in (when plotting multiple layers)}
  \item{maxnl}{integer > 0. Maximum number of layers to plot (for a multi-layer object)}
  \item{main}{character. Main plot title}
  \item{gridded}{logical. If \code{TRUE} the scatterplot is gridded (counts by cells)}
  \item{ncol}{integer. Number of columns for gridding}
  \item{nrow}{integer. Number of rows for gridding}
  \item{add}{Add points, lines, or polygons to the current plot (as an alternative to the \code{lines} or \code{points} methods}   
  \item{...}{additional graphical arguments}  
}


\examples{ 
f <- system.file("exdata/test.tif", package="terra") 
r <- rast(f)
plot(r)

d <- (r > 400) + (r > 600)
plot(d, leg.main="Title")
plot(d, leg.ext=c(178000,178200,332800,333250),
		leg.main=c("Title     \n\n"))

# character labels
x <- round(r/1000)
x <- as.factor(x)
levels(x) <- c("earth", "wind", "fire")
plot(x)


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

