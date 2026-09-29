\name{modal}

\alias{modal}
\alias{modal,SpatRaster-method}
\alias{modal,ANY-method}

\title{modal value}

\description{
Compute the mode (most frequent value). For a SpatRaster this is done for each cell across layers. For a vector (or other atomic object) it is computed from the values of \code{x} and any additional arguments in \code{...}.
}

\usage{
\S4method{modal}{SpatRaster}(x, ..., ties="first", na.rm=FALSE, filename="", overwrite=FALSE, wopt=list())

\S4method{modal}{ANY}(x, ..., ties="random", na.rm=FALSE, freq=FALSE)
}

\arguments{
  \item{x}{SpatRaster or a vector (numeric, integer, logical, character, or factor)}
  \item{...}{additional argument of the same type as \code{x}, or numeric (SpatRaster method)}
  \item{ties}{character. Indicates how to treat ties. Either "random", "lowest", "highest", "first", or "NA"}
  \item{na.rm}{logical. If \code{TRUE}, \code{NA} values are ignored. If \code{FALSE}, \code{NA} is returned if \code{x} has any \code{NA} values}
  \item{freq}{logical. If \code{TRUE}, the frequency of the modal value is returned instead of the value itself}
  \item{filename}{character. Output filename}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list with named options for writing files as in \code{\link{writeRaster}}}  
}

\value{
SpatRaster, or a single value (or its frequency if \code{freq=TRUE})
}

\examples{
r <- rast(system.file("ex/logo.tif", package="terra"))   
r <- c(r/2, r, r*2)
m <- modal(r)

modal(c(1, 2, 2, 3, 1, 2), ties="lowest")
}

\keyword{univar}
\keyword{math}
