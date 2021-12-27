\name{contour}

\docType{methods}

\alias{contour}
\alias{contour,SpatRaster-method}
\alias{as.contour}
\alias{as.contour,SpatRaster-method}

\title{Contour plot}

\description{
Contour lines of a SpatRaster. Use \code{add=TRUE} to add the lines to the current plot.  See \code{\link[graphics]{contour}} for details. 

if \code{filled=TRUE}, a new filled contour plot is made. See \code{\link[graphics]{filled.contour}} for details.

\code{as.contour} returns the contour lines as a SpatVector.
}

\usage{
\S4method{contour}{SpatRaster}(x, maxcells=100000, filled=FALSE, ...)

\S4method{as.contour}{SpatRaster}(x, maxcells=100000, ...)
}

\arguments{
  \item{x}{SpatRaster. Only the first layer is used}
  \item{maxcells}{maximum number of pixels used to create the contours}
  \item{filled}{logical. If \code{TRUE}, a \code{\link[graphics]{filled.contour}} plot is made}
  \item{...}{any argument that can be passed to \code{\link[graphics]{contour}} or  \code{\link[graphics]{filled.contour}} (graphics package)}  
}


\seealso{ 
\code{\link[terra]{plot}}
}

\examples{
r <- rast(system.file("ex/elev.tif", package="terra"))
plot(r)
contour(r, add=TRUE)

v <- as.contour(r)
plot(r)
lines(v)

contour(r, filled=TRUE, nlevels=5)

## if you want a SpatVector with contour lines
template <- disagg(rast(r), 10)
rr <- resample(r, template)
rr <- floor(rr/100) * 100
v <- as.polygons(rr)
plot(v, 1, col=terrain.colors(7))

}

\keyword{methods}
\keyword{spatial}

