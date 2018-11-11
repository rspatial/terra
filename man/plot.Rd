\name{plot}

\docType{methods}

\alias{plot}
\alias{plot,SpatRaster,missing-method}
\alias{plot,SpatLayer,missing-method}
\alias{plot,SpatRaster,numeric-method}
\alias{image}
\alias{image,SpatRaster-method}
\alias{points,SpatLayer-method}
\alias{lines,SpatLayer-method}


\title{Plot a SpatRaster or SpatLayer}

\description{
Plot (that is, make a map of) the values of a SpatRaster or SpatLayer
}

\usage{
\S4method{plot}{SpatRaster,missing}(x, y, maxpixels=100000, xlab="", ylab="", ...)
\S4method{plot}{SpatRaster,numeric}(x, y, maxpixels=100000, xlab="", ylab="", ...)
\S4method{plot}{SpatLayer,missing}(x, y, xlab="", ylab="", ...)
\S4method{image}{SpatRaster}(x, y, maxpixels=100000, xlab="", ylab="", useRaster=TRUE, ...)
\S4method{points}{SpatLayer}(x, ...)
\S4method{lines}{SpatLayer}(x, ...)
}

\arguments{
  \item{x}{SpatRaster or SpatLayer}
  \item{y}{missing or integer indicating the layer(s) to be plotted}
  \item{maxpixels}{integer > 0. Maximum number of cells to use for the plot}
  \item{xlab}{character}
  \item{ylab}{character}
  \item{useRaster}{logical; if TRUE a bitmap raster is used to plot the image instead of polygons}
  \item{...}{additional arguments for levelplot}  
}


\examples{ 
f <- system.file("exdata/test.tif", package="terra") 
r <- rast(f)
plot(r)

f <- system.file("exdata/lux.shp", package="terra")
v <- vect(f)
r <- rast(v)
values(r) <- 1:ncell(r)
image(r)
lines(v)
points(v)
}

\keyword{methods}
\keyword{spatial}

