\name{gdal}

\alias{gdal}
\alias{libVersion}
\alias{gdalCache}
\alias{clearVSIcache}
\alias{getGDALconfig}
\alias{setGDALconfig}
\alias{unloadGDALdrivers}
\alias{proj_ok}
\alias{projNetwork}


\title{GDAL version, supported file formats, cache size, and PROJ coordinate transformation control}

\description{
Set the \code{GDAL} warning level or get a \code{data.frame} with the available GDAL drivers (file formats), or, if \code{warn=NA} and \code{drivers=FALSE}, you get the version numbers of one or all of the GDAL, PROJ and GEOS libraries. 

\code{GDAL} is the software library that terra builds on to read and write spatial data and for some raster data processing. \code{PROJ} is used for transformation of coordinates ("projection") and \code{GEOS} is used for geometric operations with vector data.

The current \code{GDAL} configuration options and obtained with \code{getGDALconfig} and changed with \code{setGDALconfig}. 

\code{proj_ok} checks if the PROJ database with CRS definitions can be found.

\code{projNetwork} controls whether PROJ can access network resources for coordinate transformations. By default, network access is enabled to provide high-accuracy grid-based datum transformations where available. When disabled, PROJ falls back to \href{https://proj.org/en/stable/glossary.html#term-Ballpark-transformation}{"ballpark transformations"} (of unknown accuracy) which could lead to errors of hundreds of meters in some cases. When enabled, PROJ downloads transformation grids from \code{https://cdn.proj.org}, requiring network connectivity and valid SSL certificates. After a successful run, the files are cached locally.
}

\usage{
gdal(warn=NA, drivers=FALSE, ...)
gdalCache(size=NA)
setGDALconfig(option, value="")
getGDALconfig(option)
clearVSIcache()
libVersion(lib="all", parse=FALSE)
unloadGDALdrivers(x)
proj_ok()
projNetwork(enable, url="")
}

\arguments{
  \item{warn}{If \code{NA} and \code{drivers=FALSE}, the version of the library specified by \code{lib} is returned. Otherwise, the value should be an integer between 1 and 4 representing the level of GDAL warnings and errors that are passed to R. 1 = warnings and errors; 2 = errors only (recoverable errors as a warning); 3 = irrecoverable errors only; 4 = ignore all errors and warnings. The default setting is 2}
  \item{drivers}{logical. If \code{TRUE} a data.frame with the raster and vector data formats that are available.} 
  \item{...}{additional arguments (for backwards compatibility only)}  
  \item{size}{numeric. The new cache size in MB}  
  \item{option}{character. GDAL configuration option name, or a "name=value" string (in which case the value argument is ignored}
  \item{value}{character. value for GDAL configuration option. Use "" to reset it to its default value}
  \item{lib}{character. "gdal", "proj", or "geos", or any other value to get the versions numbers of all three}
  \item{parse}{logical. Should the version be parsed into three numerical values (major, minor and sub versions)?}
  \item{x}{character. Drivers names such as "GTiff" to be unloaded. Or "" to reload all drivers}
  \item{enable}{logical. If \code{TRUE}, enable PROJ network access for high-accuracy grid-based transformations. If \code{FALSE}, disable network access and use ballpark transformations. If missing, return the current network status}
  \item{url}{character. Optional URL for PROJ network endpoint. If empty string (default), uses PROJ's default network settings (\code{https://cdn.proj.org})}
}


\seealso{\code{\link{describe}} for file-level metadata "GDALinfo"}

\value{
character (for most functions), logical (for \code{proj_ok} and \code{projNetwork})
}

\examples{
gdal()
gdal(2)
head(gdal(drivers=TRUE))
libVersion("all", TRUE)
projNetwork()
}

\note{
While some spatial analyses may not be greatly affected by PROJ network settings (ballpark vs. grid-based transformations), the differences can be significant, especially when a transformation involves a shift in datum between different coordinate reference systems. For applications requiring high positional accuracy, ensure network access is enabled or grids are locally available. Grids can be pre-downloaded using the \code{projsync} utility or installed via system packages, such as \code{proj-data} on Ubuntu/Debian systems. Downloaded grids are cached locally and then reused for subsequent transformations. 
}

\keyword{spatial}

