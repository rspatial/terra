\name{isLonLat}

\alias{isLonLat}
\alias{isLonLat,SpatRaster-method}
\alias{isLonLat,SpatVector-method}

  
\title{Check for longitude/latitude crs}

\description{
Test whether a SpatRaster or SpatVector has a longitude/latitude coordinate reference system (CRS), or perhaps has one. That is wen the CRS is unknown (\code{""}) but the x coordinates are within -181 and 181 and the y coordinates are within -90.1 and 90.1. For a SpatRaster you can also test if it is longitude/latitude and "global" (covers all longitudes).
}

\usage{
\S4method{isLonLat}{SpatRaster}(x, perhaps=FALSE, warn=TRUE, global=FALSE, ...)

\S4method{isLonLat}{SpatVector}(x, perhaps=FALSE, warn=TRUE, ...)
}

\arguments{
\item{x}{SpatRaster or SpatVector object}
\item{perhaps}{logical. If \code{TRUE} and the crs is unknown, the method returns \code{TRUE} if the coordinates are plausible for longitude/latitude}
\item{warn}{logical. If \code{TRUE}, a warning is given if the CRS is unknown or when the CRS is longitude/latitude but the coordinates do not match that}
\item{global}{logical. If \code{TRUE}, the method tests if the raster covers all longitudes (from -180 to 180 degrees) such that the extreme columns are in fact adjacent }
\item{...}{additional arguments. None implemented}
}

\value{
logical
}

\examples{
r <- rast()
isLonLat(r)
isLonLat(r, global=TRUE)

crs(r) <- ""
isLonLat(r)
isLonLat(r, perhaps=TRUE, warn=FALSE)

crs(r) <- "+proj=lcc +lat_1=48 +lat_2=33 +lon_0=-100 +ellps=WGS84"
isLonLat(r)
}

\keyword{spatial}

