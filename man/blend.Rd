\name{blend}

\docType{methods}

\alias{blend}
\alias{blend,SpatRaster,SpatRaster-method}
\alias{blend,SpatRasterCollection,missing-method}

\title{ 
Distance-weighted raster blending
}

\description{
Combine rasters with smooth gradients in overlappoing zones by weighting each raster's contribution by the distance from the cell to the raster's nearest edge. The result is independent of the order of the rasters.

The SpatRasters must have the same origin and spatial resolution. 
}

\usage{
\S4method{blend}{SpatRaster,SpatRaster}(x, y, ..., filename="", overwrite=FALSE, wopt=list())

\S4method{blend}{SpatRasterCollection,missing}(x, filename="", ...)
}

\arguments{
  \item{x}{SpatRaster or SpatRasterCollection}
  \item{y}{SpatRaster if \code{x} is a SpatRaster, otherwise missing}
  \item{...}{additional SpatRasters if \code{x} is a SpatRaster, otherwise options for writing files as in \code{\link{writeRaster}}}
  \item{filename}{character. Output filename}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list with named options for writing files as in \code{\link{writeRaster}}}
} 


\value{
SpatRaster
}

\seealso{
\code{\link{mosaic}, \link{merge}}
}


\examples{

x <- rast(xmin=-110, xmax=-60, ymin=40, ymax=70, res=1, vals=1)
y <- rast(xmin=-95, xmax=-45, ymax=60, ymin=30, res=1, vals=2)
z <- rast(xmin=-80, xmax=-30, ymax=50, ymin=20, res=1, vals=3)

m1 <- blend(x, y, z)

# with many SpatRasters, make a SpatRasterCollection from a list
rlist <- list(x, y, z)
rsrc <- sprc(rlist)

m2 <- blend(rsrc)
}

\keyword{methods}
\keyword{spatial}

