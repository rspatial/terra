\name{scatterplot}

\docType{methods}

\alias{plot,SpatRaster,SpatRaster-method}


\title{Scatterplot of two SpatRaster layers}

\description{
Scatterplot of the values of two SpatRaster layers
}

\usage{
\S4method{plot}{SpatRaster,SpatRaster}(x, y, maxcell=100000, warn=TRUE, nc, nr, 
   maxnl=16, smooth=FALSE, gridded=FALSE, ncol=25, nrow=25, ...)

}

\arguments{
  \item{x}{SpatRaster}
  \item{y}{SpatRaster}
  \item{maxcell}{positive integer. Maximum number of cells to use for the plot}
  \item{nc}{positive integer. Optional. The number of columns to divide the plotting device in (when plotting multiple layers)}
  \item{nr}{positive integer. Optional. The number of rows to divide the plotting device in (when plotting multiple layers)}
  \item{maxnl}{positive integer. Maximum number of layers to plot (for multi-layer objects)}
  \item{smooth}{logical. If \code{TRUE} show a smooth scatterplot}
  \item{gridded}{logical. If \code{TRUE} the scatterplot is gridded (counts by cells)}
  \item{warn}{boolean. Show a warning if a sample of the pixels is used (for scatterplot only)}

  \item{ncol}{positive integer. Number of columns for gridding}
  \item{nrow}{positive integer. Number of rows for gridding}
  \item{...}{additional graphical arguments}  
}


\examples{ 
s <- rast(system.file("ex/logo.tif", package="terra"))   
plot(s[[1]], s[[2]])
plot(s, sqrt(s[[3:1]]))
}

\keyword{methods}
\keyword{spatial}

