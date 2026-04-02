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

With the "blend" function, smooth gradients in overlappoing zones are created by weighting each raster's contribution by the distance from the cell to the raster's nearest edge. 

The SpatRasters must have the same origin and spatial resolution. 

If rasters do not overlap, or if their values do not need to be computed, you can use \code{\link{merge}} instead. 
}

\usage{
\S4method{mosaic}{SpatRaster,SpatRaster}(x, y, ..., fun="mean", resample=FALSE, method="", filename="", overwrite=FALSE, wopt=list())

\S4method{mosaic}{SpatRasterCollection,missing}(x, fun="mean", resample=FALSE, method="", filename="", ...)
}

\arguments{
  \item{x}{SpatRaster or SpatRasterCollection}
  \item{y}{SpatRaster if \code{x} is a SpatRaster, otherwise missing}
  \item{...}{additional SpatRasters}
  \item{fun}{character. One of "blend", "mean", "median", "min", "max", "modal", "sum", "first", "last"}
  \item{resample}{logical. If \code{TRUE} input rasters are resampled if they do not align (same origin and resolution) with the first raster}
  \item{method}{character. The resampling method used. One of "nearest", "bilinear", "cubic", "cubicspline", "lanczos", "average", "mode" as in \code{\link{resample}}. If \code{NULL}, "nearest" is used for categorical rasters and "bilinear" for other rasters}   
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

# with a SpatRasterCollection
spc <- sprc(list(x, y, z))

m <- mosaic(spc)

b <- mosaic(spc, fun="blend")
}

\keyword{methods}
\keyword{spatial}

