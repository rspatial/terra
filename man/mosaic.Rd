\name{mosaic}

\docType{methods}

\alias{mosaic}
\alias{mosaic,SpatRaster,SpatRaster-method}
\alias{mosaic,SpatRasterCollection,missing-method}

\title{ 
mosaic SpatRasters
}

\description{
Combine adjacent and (partly) overlapping SpatRasters to form a single new SpatRaster. Values in overlapping cells are averaged (by default) or can be computed with another function.

The SpatRasters must have the same origin and spatial resolution. 

This method is similar to the simpler, but much faster, \code{\link{merge}} method. 
}

\usage{
\S4method{mosaic}{SpatRaster,SpatRaster}(x, y, ..., fun="mean", filename="", overwrite=FALSE, wopt=list())

\S4method{mosaic}{SpatRasterCollection,missing}(x, fun="mean", filename="", ...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{y}{object of same class as \code{x}}
  \item{...}{additional SpatRasters}
  \item{fun}{character. One of "mean", "median", "min", "max", "modal", "sum", "first", "last"}
  \item{filename}{character. Output filename}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list with named options for writing files as in \code{\link{writeRaster}}}
} 


\value{
SpatRaster
}

\seealso{
\code{\link{merge}}
}


\examples{

x <- rast(xmin=-110, xmax=-60, ymin=40, ymax=70, res=1, vals=1)
y <- rast(xmin=-95, xmax=-45, ymax=60, ymin=30, res=1, vals=2)
z <- rast(xmin=-80, xmax=-30, ymax=50, ymin=20, res=1, vals=3)

m1 <- mosaic(x, y, z)

m2 <- mosaic(z, y, x)

# with many SpatRasters, make a SpatRasterCollection from a list
rlist <- list(x, y, z)
rsrc <- sprc(rlist)

m <- mosaic(rsrc)
}

\keyword{methods}
\keyword{spatial}

