\name{plotVector}

\docType{methods}

\alias{plot,SpatVector,missing-method}
\alias{plot,SpatVector,character-method}
\alias{plot,SpatExtent,missing-method}
\alias{lines,SpatVector-method}
\alias{lines,SpatExtent-method}
\alias{lines,SpatRaster-method}
\alias{polys}
\alias{polys,SpatVector-method}
\alias{points,SpatVector-method}
\alias{points,SpatExtent-method}


\title{Plot a SpatVector}

\description{
Plot SpatVector data to make a map
}

\usage{
\S4method{plot}{SpatVector,missing}(x, y, col, 
   axes=TRUE, add=FALSE, ...)

\S4method{plot}{SpatVector,character}(x, y, col, type,
   mar=c(5.1, 4.1, 4.1, 7.1), axes=TRUE, legend=TRUE, add=FALSE, ...)

\S4method{plot}{SpatExtent,missing}(x, y, col, ...)

\S4method{points}{SpatVector}(x, col, ...)

\S4method{lines}{SpatVector}(x, col, ...)

\S4method{polys}{SpatVector}(x, col, border="black", ...)

\S4method{points}{SpatExtent}(x, col, ...)

\S4method{lines}{SpatExtent}(x, col, ...)

}

\arguments{
  \item{x}{SpatVector}
  \item{y}{missing or positive integer indicating the layer(s) to be plotted, or the name of the layer to be mapped}
  \item{col}{character. Colors}
  \item{type}{character. Type of map/legend. One of "continuous", "classes", or "interval"} 
  \item{axes}{logical. Draw axes?}
  \item{legend}{logical. Draw a legend?}
  \item{mar}{numeric vector of lenght 4 to set the margins of the plot (to make space for the legend)}
  \item{border}{character. color(s) of the polygon borders. Use \code{NULL} or \code{NA} to not draw a border}
  \item{add}{logical. If \code{TRUE} add the object to the current plot (for points, lines, and polygons). This allows you to add vector data to a map, specifying a field, and requesting a legend.}
  \item{...}{additional graphical arguments such as \code{lwd}, \code{cex} and \code{pch}}  
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

