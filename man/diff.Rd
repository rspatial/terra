\name{diff}

\alias{diff}
\alias{diff,SpatRaster-method}
  
\title{Lagged differences}

\description{
Compute the difference between consecutive layers in a SpatRaster.
}

\usage{
\S4method{diff}{SpatRaster}(x, lag=1, filename="", ...)
}


\arguments{
  \item{x}{SpatRaster}
  \item{lag}{positive integer indicating which lag to use}
  \item{filename}{character. Output filename}
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
}

\value{
SpatRaster
}


\examples{
s <- rast(system.file("ex/logo.tif", package="terra"))   
d <- diff(s)
}

\keyword{methods}
\keyword{spatial}
