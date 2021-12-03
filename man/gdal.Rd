\name{gdal}

\alias{gdal}
\alias{gdalCache}

\title{GDAL version, supported file formats, and cache size}

\description{
Set the \code{GDAL} warning level or get a \code{data.frame} with the available GDAL drivers (file formats), or get the GDAL version (if \code{warn=NA} and \code{drivers=FALSE}, the default. 

\code{GDAL} is the software library that terra builds on to read and write spatial data and for some raster data processing.
}

\usage{
gdal(warn=NA, drivers=FALSE)
gdalCache(size=NA)
}

\arguments{
  \item{warn}{ignored if \code{NA}. Otherwise, the value should be an integer between 1 and 4 representing the level of GDAL warnings and errors that are passed to R. 1 = warnings and errors; 2 = errors only (recoverable errors as a warning); 3 = irrecoverable errors only; 4 = ignore all errors and warnings. The default setting is 3}
  \item{drivers}{logical. If \code{TRUE} a data.frame with the raster and vector data formats that are available.}  
  \item{size}{numeric. The new cache size in MB}  
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

