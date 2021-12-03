\name{is.bool}

\docType{methods}

\alias{is.bool}
\alias{is.bool,SpatRaster-method}
\alias{as.bool}
\alias{as.bool,SpatRaster-method}
\alias{isTRUE,SpatRaster-method}
\alias{isFALSE,SpatRaster-method}
\alias{as.logical,SpatRaster-method}
\alias{is.int}
\alias{is.int,SpatRaster-method}
\alias{as.int}
\alias{as.int,SpatRaster-method}
\alias{as.integer,SpatRaster-method}

  
\title{Raster value types}

\description{
The values in a SpatRaster layer are by default numeric, but they can also be logical (Boolean), integer, or categorical

Note that \code{as.bool} and \code{as.int} return a new \code{SpatRaster}, whereas \code{is.bool} and \code{is.int} return a \code{logical} value for each layer. For a \code{SpatRaster}, \code{isTRUE} is equivalent to \code{as.bool}, \code{isFALSE} is equivalent to \code{!as.bool}, code{as.integer} is the same as \code{as.int} and \code{as.logical} is the same as \code{as.bool}
}

\usage{
\S4method{is.bool}{SpatRaster}(x)
\S4method{as.bool}{SpatRaster}(x, filename, ...)
\S4method{is.int}{SpatRaster}(x)
\S4method{as.int}{SpatRaster}(x, filename, ...)
}


\arguments{
  \item{x}{SpatRaster}
  \item{filename}{character. Output filename}
  \item{...}{list with named options for writing files as in \code{\link{writeRaster}}}
}  


\value{
logical or SpatRaster
}



\examples{
r <- rast(nrows=10, ncols=10, vals=1:100)
is.bool(r)
z <- as.bool(r)
is.bool(z)

x <- r > 25
is.bool(x)

rr <- r/2
is.int(rr)
is.int(round(rr))
}

\keyword{methods}
\keyword{spatial}

