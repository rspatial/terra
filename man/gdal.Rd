\name{gdal}

\alias{gdal}
\alias{gdalCache}
\alias{getGDALconfig}
\alias{setGDALconfig}

\title{GDAL version, supported file formats, and cache size}

\description{
Set the \code{GDAL} warning level or get a \code{data.frame} with the available GDAL drivers (file formats), or, if \code{warn=NA} and \code{drivers=FALSE}, you get the version numbers of one or all of the GDAL, PROJ and GEOS libraries. 

\code{GDAL} is the software library that terra builds on to read and write spatial data and for some raster data processing. \code{PROJ} is used for transformation of coordinates ("projection") and \code{GEOS} is used for geometric operations with vector data.
}

\usage{
gdal(warn=NA, drivers=FALSE, lib="gdal")
gdalCache(size=NA)
setGDALconfig(option, value="")
getGDALconfig(option)
}

\arguments{
  \item{warn}{If \code{NA} and \code{drivers=FALSE}, the version of the library specified by \code{lib} is returned. Otherwise, the value should be an integer between 1 and 4 representing the level of GDAL warnings and errors that are passed to R. 1 = warnings and errors; 2 = errors only (recoverable errors as a warning); 3 = irrecoverable errors only; 4 = ignore all errors and warnings. The default setting is 3}
  \item{drivers}{logical. If \code{TRUE} a data.frame with the raster and vector data formats that are available.} 
  \item{lib}{character. "gdal", "proj", or "geos", or any other value to get the versions numbers of all three}
  \item{size}{numeric. The new cache size in MB}  
  \item{option}{character. GDAL configuration option name, or a "name=value" string (in which case the value argument is ignored}
  \item{value}{character. value for GDAL configuration option. Use "" to reset it to its default value}
}


\seealso{\code{\link{describe}} for file-level metadata "GDALinfo"}

\value{
character
}

\examples{
gdal()
gdal(2)
head(gdal(drivers=TRUE))
}

\keyword{spatial}

