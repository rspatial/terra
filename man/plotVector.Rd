\name{plotVector}

\docType{methods}

\alias{plot,SpatVector,missing-method}
\alias{plot,SpatVector,character-method}
\alias{plot,SpatExtent,missing-method}
\alias{lines,SpatVector-method}
\alias{lines,SpatExtent-method}
\alias{lines,SpatRaster-method}
\alias{points,SpatVector-method}
\alias{points,SpatExtent-method}


\title{Plot a SpatVector}

\description{
Plot (that is, make a map of) using a SpatVector
}

\usage{
\S4method{plot}{SpatVector,missing}(x, y, col, 
   axes=TRUE, add=FALSE, ...)

\S4method{plot}{SpatVector,character}(x, y, col, type,
   mar=c(5.1, 4.1, 4.1, 7.1), axes=TRUE, add=FALSE, ...)

\S4method{plot}{SpatExtent,missing}(x, y, col, ...)

\S4method{points}{SpatVector}(x, col, ...)

\S4method{lines}{SpatVector}(x, col, ...)

\S4method{points}{SpatExtent}(x, col, ...)

\S4method{lines}{SpatExtent}(x, col, ...)

}

\arguments{
  \item{x}{SpatVector}
  \item{y}{missing or positive integer indicating the layer(s) to be plotted, or the name of the layer to be mapped}
  \item{col}{character. Colors}
  \item{type}{character. Type of map/legend. One of "continuous", "classes", or "internval"}
  
  \item{axes}{logical. Draw axes?}
  
  \item{mar}{numeric vector of lenght 4 to set the margins of the plot (to make space for the legend)}
  \item{add}{logical. If \code{TRUE} add the object to the current plot (for points, lines, and polygons  this is an alternative to the \code{lines} or \code{points} methods}   
  \item{...}{additional graphical arguments}  
}


\examples{ 
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
plot(v)
plot(v, "ID_1")

r <- rast(v)
values(r) <- 1:ncell(r)
plot(r)
lines(v)
points(v)
}

\keyword{methods}
\keyword{spatial}

