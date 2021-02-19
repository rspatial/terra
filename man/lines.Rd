\name{lines}

\docType{methods}

\alias{points}
\alias{lines}
\alias{polys}

\alias{lines,SpatVector-method}
\alias{lines,SpatRaster-method}
\alias{polys,SpatVector-method}
\alias{points,SpatVector-method}

\alias{lines,SpatExtent-method}
\alias{points,SpatExtent-method}


\title{Add SpatVector data to a map}

\description{
Add SpatVector data to a plot (map) with \code{points}, \code{lines}, or \code{polys}. 

These are simpler alternatives for \code{\link[terra:plot]{plot(x, add=TRUE)}}
}

\usage{
\S4method{points}{SpatVector}(x, col, cex=1, pch=20, ...)

\S4method{lines}{SpatVector}(x, y=NULL, col, lwd=1, lty=1, arrows=FALSE, ...)

\S4method{polys}{SpatVector}(x, col, border="black", lwd=1, lty=1, ...)

\S4method{points}{SpatExtent}(x, col, ...)

\S4method{lines}{SpatExtent}(x, col, ...)
}

\arguments{
  \item{x}{SpatVector or SpatExtent}
  \item{y}{missing or SpatVector. If both \code{x} and \code{y} have point geometry and the same number of rows, lines are drawn between pairs of points}
  \item{col}{character. Colors}
  \item{border}{character. color(s) of the polygon borders. Use \code{NULL} or \code{NA} to not draw a border}
  \item{cex}{numeric. point size magnifier. See \code{\link{par}}}
  \item{pch}{positive integer, line type. See \code{\link{points}}}
  \item{lwd}{numeric, line-width. See \code{\link{par}}}
  \item{lty}{positive integer, line type. See \code{\link{par}}}
  \item{arrows}{logical. If \code{TRUE} and \code{y} is a SpatVector, arrows are drawn intead of lines. See \code{?\link{arrows}} for additional arguments}
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

