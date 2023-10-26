\name{is.lonlat}

\alias{is.lonlat}
\alias{is.lonlat,SpatRaster-method}
\alias{is.lonlat,SpatVector-method}
\alias{is.lonlat,character-method}

  
\title{Check for longitude/latitude crs}

\description{
Test whether a SpatRaster or SpatVector has a longitude/latitude coordinate reference system (CRS), or perhaps has one. That is, when the CRS is unknown (\code{""}) but the x coordinates are within -181 and 181 and the y coordinates are within -90.1 and 90.1. For a SpatRaster you can also test if it has a longitude/latitude CRS and it is "global" (covers all longitudes).

A warning is given if the CRS is missing or if it is specified as longitude/latitude but the coordinates do not match that.
}

\usage{
\S4method{is.lonlat}{SpatRaster}(x, perhaps=FALSE, warn=TRUE, global=FALSE)

\S4method{is.lonlat}{SpatVector}(x, perhaps=FALSE, warn=TRUE)

\S4method{is.lonlat}{character}(x, perhaps=FALSE, warn=TRUE)
}

\arguments{
\item{x}{SpatRaster or SpatVector}
\item{perhaps}{logical. If \code{TRUE} and the CRS is unknown, the method returns \code{TRUE} if the coordinates are plausible for longitude/latitude}
\item{warn}{logical. If \code{TRUE}, a warning is given if the CRS is unknown but assumed to be lon/lat and  \code{perhaps=TRUE}}
\item{global}{logical. If \code{TRUE}, the method tests if the raster covers all longitudes (from -180 to 180 degrees) such that the extreme columns are in fact adjacent}
}

\value{
logical or NA
}

\examples{
r <- rast()
is.lonlat(r)
is.lonlat(r, global=TRUE)

crs(r) <- ""
is.lonlat(r)
is.lonlat(r, perhaps=TRUE, warn=FALSE)

crs(r) <- "+proj=lcc +lat_1=48 +lat_2=33 +lon_0=-100 +ellps=WGS84"
is.lonlat(r)
}

\keyword{spatial}

