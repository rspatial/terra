\name{isLonLat}

\alias{isLonLat}
\alias{isLonLat,SpatRaster-method}
\alias{isLonLat,SpatVector-method}
\alias{isGlobalLonLat}
\alias{isGlobalLonLat,SpatRaster-method}
\alias{couldBeLonLat}
\alias{couldBeLonLat,SpatRaster-method}
\alias{couldBeLonLat,SpatVector-method}

  
\title{Check for longitude/latitude crs}

\description{
Test whether a SpatRaster or SpatVector has a longitude/latitude coordinate reference system (CRS). \code{couldBeLonLat} returns \code{TRUE} if the CRS is unknown (\code{""}) but the x coordinates are within -181 and 181 and the y coordinates are within -90.1 and 90.1.
}

\usage{
\S4method{isLonLat}{SpatRaster}(x, ...)
\S4method{isLonLat}{SpatVector}(x, ...)
\S4method{isGlobalLonLat}{SpatRaster}(x, ...)
\S4method{couldBeLonLat}{SpatRaster}(x, warn=TRUE, ...)
\S4method{couldBeLonLat}{SpatVector}(x, warn=TRUE, ...)
}

\arguments{
\item{x}{SpatRaster or SpatVector object}
\item{warn}{logical. If \code{TRUE}, a warning is given if the CRS is unknown or when the CRS is longitude/latitude but the coordinates do not match that}
\item{...}{additional arguments. None implemented}
}

\value{
Logical
}

\examples{
r <- rast()
isLonLat(r)
crs(r) <- "+proj=lcc +lat_1=48 +lat_2=33 +lon_0=-100 +ellps=WGS84"
isLonLat(r)
}

\keyword{spatial}

