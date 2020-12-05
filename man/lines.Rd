\name{lines}

\docType{methods}

\alias{lines,SpatVector-method}
\alias{lines,SpatRaster-method}
\alias{polys}
\alias{polys,SpatVector-method}
\alias{points,SpatVector-method}

\alias{lines,SpatExtent-method}
\alias{points,SpatExtent-method}


\title{Add SpatVector data to a map}

\description{
Add SpatVector data to a plot (map) with \code{points}, \code{lines}, or \code{polys}. 

These are simpler alternatives for \code{plot(x, add=TRUE)}
}

\usage{
\S4method{points}{SpatVector}(x, col, ...)

\S4method{lines}{SpatVector}(x, col, ...)

\S4method{polys}{SpatVector}(x, col, border="black", ...)

\S4method{points}{SpatExtent}(x, col, ...)

\S4method{lines}{SpatExtent}(x, col, ...)
}

\arguments{
  \item{x}{SpatVector or SpatExtent}
  \item{col}{character. Colors}
  \item{border}{character. color(s) of the polygon borders. Use \code{NULL} or \code{NA} to not draw a border}
  \item{...}{additional graphical arguments such as \code{lwd}, \code{cex} and \code{pch}}  
}


\examples{ 
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

