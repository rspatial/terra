\name{image}

\docType{methods}

\alias{image}
\alias{image,SpatRaster-method}

\title{SpatRaster image method}

\description{
Plot (make a map of) the values of a SpatRaster via \code{\link[graphics]{image}}. See \code{\link[terra]{plot}} if you need more fancy options such as a legend.
}

\usage{
\S4method{image}{SpatRaster}(x, y=1, maxcell=500000, ...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{y}{positive integer indicating the layer to be plotted, or a character indicating the name of the layer}
  \item{maxcell}{positive integer. Maximum number of cells to use for the plot}
  \item{...}{additional arguments as for \code{graphics::\link[graphics]{image}} }
}


\seealso{
 \code{\link[terra]{plot}}  
}

\examples{ 
f <- system.file("ex/elev.tif", package="terra") 
r <- rast(f)
image(r)
image(r, col=rainbow(24))
}

\keyword{methods}
\keyword{spatial}

