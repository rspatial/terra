\name{gdal}

\alias{gdal}
\alias{gdal}

\title{gdal utilities}

\description{
Set the GDAL warning level or get a data.frame with the available GDAL drivers (file formats), or get the GDAL version (if \code{warn=NA} and \code{drivers=FALSE}, the default. 
}

\usage{
gdal(warn=NA, drivers=FALSE)
}

\arguments{
  \item{warn}{ignored if \code{NA}. Otherwise, the value should be an integer between 1 and 4 representing the level of GDAL warnings and errors that are passed to R. 1 = warnings and errors; 2 = errors only (recoverable errors as a warning); 3 = irrecoverable errors only; 4 = ignore all errors and warnings. The default setting is 3}
  \item{drivers}{logical. If \code{TRUE} a data.frame with the raster and vector data formats that are available.}  
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

