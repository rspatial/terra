\name{compareGeom}

\alias{compareGeom}
\alias{compareGeom,SpatRaster,SpatRaster-method}
\alias{compareGeom,SpatRaster,SpatRasterCollection-method}

\alias{compareGeom,SpatRaster,list-method}
\alias{compareGeom,SpatRasterCollection,missing-method}

\alias{compareGeom,SpatVector,SpatVector-method}
\alias{compareGeom,SpatVector,missing-method}
  
\title{Compare geometries of SpatRasters}

\description{
Evaluate whether two SpatRasters have the same extent, number of rows and columns, projection, resolution, and origin (or a subset of these comparisons). 

Or evaluate whether two SpatVectors have the same geometries, or whether a SpatVector has duplicated geometries. 
}

\usage{
\S4method{compareGeom}{SpatRaster,SpatRaster}(x, y, ..., lyrs=FALSE, crs=TRUE, warncrs=FALSE, ext=TRUE,
	rowcol=TRUE, res=FALSE, stopOnError=TRUE, messages=FALSE)

\S4method{compareGeom}{SpatVector,SpatVector}(x, y, tolerance=0)

\S4method{compareGeom}{SpatVector,missing}(x, y, tolerance=0)
}

\arguments{
  \item{x}{SpatRaster or SpatVector}
  \item{y}{Same as \code{x}. If \code{x} is a SpatRaster, \code{y} can also be a list of SpatRasters. If \code{x} is a SpatVector, \code{y} can be missing}
  \item{...}{Additional SpatRasters}
  \item{lyrs}{logical. If \code{TRUE}, the number of layers is compared}
  \item{crs}{logical. If \code{TRUE}, coordinate reference systems are compared}
  \item{warncrs}{logical. If \code{TRUE}, a warning is given if the crs is different (instead of an error)}
  \item{ext}{logical. If \code{TRUE}, bounding boxes are compared}
  \item{rowcol}{logical. If \code{TRUE}, number of rows and columns of the objects are compared}
  \item{res}{logical. If \code{TRUE}, resolutions are compared (redundant when checking extent and rowcol)}
  \item{stopOnError}{logical. If \code{TRUE}, code execution stops if raster do not match}
  \item{messages}{logical. If \code{TRUE}, warning/error messages are printed even if \code{stopOnError=FALSE}}
  \item{tolerance}{numeric}
}

\value{
logical (SpatRaster) or matrix of logical (SpatVector)
} 

\examples{
r1 <- rast()
r2 <- rast()
r3 <- rast()
compareGeom(r1, r2, r3)
nrow(r3) <- 10


\dontrun{
compareGeom(r1, r3)
}
}

\keyword{spatial}

