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
\alias{is.factor}
\alias{is.factor,SpatRaster-method}
\alias{as.factor}
\alias{as.factor,SpatRaster-method}

  
\title{Raster value types}

\description{
The values in a SpatRaster layer are by default numeric, but they can also be set to be logical (Boolean), integer, or categorical (factor). 

For a \code{SpatRaster}, \code{as.logical} and \code{isTRUE} is equivalent to \code{as.bool}. \code{isFALSE} is equivalent to \code{!as.bool}, and \code{as.integer} is the same as \code{as.int}.

\code{as.bool} and \code{as.int} force the values into the correct range (e.g. whole integers) but in-memory cell values are still stored as numeric. They will behave like the assigned types, though, and will be written to files with that data type (if the file type supports it). 

See \code{\link{levels}} and \code{\link{cats}} to create categorical layers by setting labels. 
}

\usage{
\S4method{is.bool}{SpatRaster}(x)

\S4method{as.bool}{SpatRaster}(x, filename, ...)

\S4method{is.int}{SpatRaster}(x)

\S4method{as.int}{SpatRaster}(x, filename, ...)

\S4method{is.factor}{SpatRaster}(x)

\S4method{as.factor}{SpatRaster}(x)
}


\arguments{
  \item{x}{SpatRaster}
  \item{filename}{character. Output filename}
  \item{...}{list with named options for writing files as in \code{\link{writeRaster}}}
}  

\seealso{
\code{\link{levels}} and \code{\link{cats}} to create categorical layers (and set labels). 
}


\value{
The \code{as.*} methods return a new \code{SpatRaster}, whereas the \code{is.*} methods return a \code{logical} value for each layer in \code{x}. 
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

