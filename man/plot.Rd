\name{plot}

\docType{methods}

\alias{plot}
\alias{plot,SpatRaster,missing-method}
\alias{plot,SpatRaster,numeric-method}
\alias{image}
\alias{image,SpatRaster-method}


\title{Plot a SpatRaster object}

\description{
Plot (that is, make a map of) the values of a SpatRaster object
}

\usage{
\S4method{plot}{SpatRaster,missing}(x, y, maxpixels=500000, xlab="", ylab="", ...)
\S4method{image}{SpatRaster}(x, maxpixels=500000, xlab="", ylab="", useRaster=TRUE, ...)
}

\arguments{
  \item{x}{Raster* object}
  \item{y}{missing or integer indicating the layer(s) to be plotted}
  \item{maxpixels}{integer > 0. Maximum number of cells to use for the plot}
  \item{xlab}{character}
  \item{ylab}{character}
  \item{useRaster}{logical; if TRUE a bitmap raster is used to plot the image instead of polygons}
  \item{...}{additional arguments for levelplot}  
}


\examples{ 
f <- system.file("external/test.grd", package="terra") 
r <- rast(f)
plot(r)
}

\keyword{methods}
\keyword{spatial}

