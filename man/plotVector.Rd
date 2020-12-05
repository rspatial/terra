\name{plotVector}

\docType{methods}

\alias{plot,SpatVector,missing-method}
\alias{plot,SpatVector,character-method}
\alias{plot,SpatVector,numeric-method}
\alias{plot,SpatExtent,missing-method}


\title{Plot a SpatVector}

\description{
Plot SpatVector data to make a map
}

\usage{
\S4method{plot}{SpatVector,missing}(x, y, col, 
   axes=TRUE, add=FALSE, pax=list(), ...)

\S4method{plot}{SpatVector,character}(x, y, col, type,
   mar=c(5.1, 4.1, 4.1, 7.1), axes=TRUE, legend=TRUE, add=FALSE, 
   plg=list(), pax=list(), ...)

\S4method{plot}{SpatVector,numeric}(x, y, col, type,
   mar=c(5.1, 4.1, 4.1, 7.1), axes=TRUE, legend=TRUE, add=FALSE, 
   plg=list(), pax=list(), ...)

\S4method{plot}{SpatExtent,missing}(x, col, ...)
}

\arguments{
  \item{x}{SpatVector}
  \item{y}{missing or positive integer indicating the layer(s) to be plotted, or the name of the layer}
  \item{col}{character. Colors}
  \item{type}{character. Type of map/legend. One of "continuous", "classes", or "interval"} 
  \item{axes}{logical. Draw axes?}
  \item{legend}{logical. Draw a legend?}
  \item{mar}{numeric vector of lenght 4 to set the margins of the plot (to make space for the legend)}
  \item{add}{logical. If \code{TRUE} add the object to the current plot}
  \item{plg}{list. Parameters for the legend}
  \item{pax}{list. Parameters for drawing axes}  
  \item{...}{additional graphical arguments such as \code{lwd}, \code{cex} and \code{pch}}  
}


\examples{ 
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
plot(v)
plot(v, density=1:12, angle=seq(18, 360, 20))

plot(v, "NAME_2", col=rainbow(12), border=c("gray", "blue"), lwd=3, type="classes")
}

\keyword{methods}
\keyword{spatial}

