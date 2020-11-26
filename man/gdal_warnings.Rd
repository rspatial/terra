\name{gdal}

\alias{gdal_version}
\alias{gdal_warnings}

\title{gdal utilities}

\description{
Get the GDAL version or set the GDAL warning level. 
}


\usage{
gdal_version()

gdal_warnings(level)
}

\arguments{
  \item{level}{integer between 1 and 4 representing the level of GDAL warnings and errors that are passed to R. 1 = warnings and errors; 2 = errors only (recoverable errors as a warning); 3 = irrecoverable errors only; 4 = ignore all errors and warnings. The default setting is 3}
}


\value{
character (gdal_version)
}


\examples{
gdal_version()
}

\keyword{spatial}

