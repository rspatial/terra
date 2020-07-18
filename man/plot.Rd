\name{plot}

\docType{methods}

\alias{plot}
\alias{plot,SpatRaster,missing-method}
\alias{plot,SpatRaster,numeric-method}
\alias{plot,SpatRaster,SpatRaster-method}

\alias{plot,SpatVector,missing-method}
\alias{plot,SpatVector,character-method}
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
\S4method{plot}{SpatRaster,numeric}(x, y, col, maxcell=100000, mar=c(5.1, 4.1, 4.1, 6.6), 
   leg.levels=5, leg.shrink=c(0,0), leg.main=NULL, leg.main.cex=1, leg.ext=NULL,
   digits, useRaster=TRUE, zlim, xlab="", ylab="", axes=TRUE, add=FALSE, ...)

\S4method{plot}{SpatRaster,missing}(x, y, maxcell=100000, nc, nr, main, 
   maxnl=16, add=FALSE, ...)

\S4method{plot}{SpatRaster,SpatRaster}(x, y, maxcell=100000, warn=TRUE, nc, nr, 
   maxnl=16, gridded=FALSE, ncol=25, nrow=25, ...)

\S4method{plot}{SpatVector,missing}(x, y, col=NULL, xlab="", ylab="", 
   axes=TRUE, add=FALSE, border="black", ...)

\S4method{plot}{SpatVector,character}(x, y, col=topo.colors(100), xlab="", ylab="", 
   axes=TRUE, add=FALSE, mar=c(5.1, 4.1, 4.1, 7.1), leg.ext=NULL, leg.type=NULL, leg.levels=5, digits, ...)

\S4method{plot}{SpatExtent,missing}(x, y, col=NULL, xlab="", ylab="", axes=TRUE, ...)

\S4method{image}{SpatRaster}(x, y, maxcell=100000, xlab="", ylab="", ...)

\S4method{points}{SpatVector}(x, col=NULL, ...)

\S4method{lines}{SpatVector}(x, col=NULL, ...)

\S4method{points}{SpatExtent}(x, col=NULL, ...)

\S4method{lines}{SpatExtent}(x, col=NULL, ...)
}

\arguments{
  \item{x}{SpatRaster or SpatVector}
  \item{y}{missing or positive integer indicating the layer(s) to be plotted, or the name of the layer to be mapped}
  \item{col}{character. Colors}
  \item{maxcell}{positive integer. Maximum number of cells to use for the plot}
  \item{mar}{numeric vector of lenght 4 to set the margins of the plot}
  \item{leg.levels}{positive integer. Number of labels on the legend}
  \item{leg.shrink}{numeric. Relative shrinkage of legend from top and bottom}
  \item{leg.main}{character. Legend title}
  \item{leg.main.cex}{character. Legend title size}
  \item{leg.ext}{SpatExtent or numeric vector indicating the area where to plot the legend. The vector should have a lenght of 4 (xmin, xmax, ymin, ymax)}
  \item{leg.type}{character. Legend type}

  \item{digits}{integer >= 0}
  \item{zlim}{range of allowed values}
  \item{xlab}{character}
  \item{ylab}{character}
  \item{axes}{logical}
  \item{useRaster}{logical; If \code{TRUE} plotting is faster because a bitmap raster is used to plot the image instead of polygons. See \code{\link[graphics]{image}}}
  \item{nc}{positive integer. Optional. The number of columns to divide the plotting device in (when plotting multiple layers)}
  \item{nr}{positive integer. Optional. The number of rows to divide the plotting device in (when plotting multiple layers)}
  \item{maxnl}{postive integer. Maximum number of layers to plot (for a multi-layer object)}
  \item{main}{character. Main plot title}
  \item{gridded}{logical. If \code{TRUE} the scatterplot is gridded (counts by cells)}
  \item{warn}{boolean. Show a warning if a sample of the pixels is used (for scatterplot only)}

  \item{ncol}{positive integer. Number of columns for gridding}
  \item{nrow}{positive integer. Number of rows for gridding}
  \item{add}{logical. If \code{TRUE} add the object to the current plot (for points, lines, and polygons  this is an alternative to the \code{lines} or \code{points} methods}   
  \item{border}{character to select the border (polygon outline) color}
  \item{...}{additional graphical arguments}  
}


\examples{ 
f <- system.file("ex/test.tif", package="terra") 
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


f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
r <- rast(v)
values(r) <- 1:ncell(r)
plot(r)
lines(v)
points(v)
}

\keyword{methods}
\keyword{spatial}

