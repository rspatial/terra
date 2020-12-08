\name{gdal}

\alias{gdal}
\alias{gdal}

\title{gdal utilities}

\description{
Get the GDAL version or set the GDAL warning level. 
}

\usage{
gdal(warn=NA)
}

\arguments{
  \item{warn}{if \code{NA} the GDAL version used is returned. Otherwise, the value should be an integer between 1 and 4 representing the level of GDAL warnings and errors that are passed to R. 1 = warnings and errors; 2 = errors only (recoverable errors as a warning); 3 = irrecoverable errors only; 4 = ignore all errors and warnings. The default setting is 3}
}


\seealso{\code{\link{describe}} for GDALinfo output}

\value{
character
}


\examples{
gdal()
gdal(3)
}

\keyword{spatial}

