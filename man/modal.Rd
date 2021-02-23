\name{modal}

\alias{modal}
\alias{modal,SpatRaster-method}

\title{modal value}

\description{
Compute the mode for each cell across the layers of a SpatRaster. The mode, or modal value, is the most frequent value in a set of values. 
}

\usage{
\S4method{modal}{SpatRaster}(x, ..., ties="first", na.rm=FALSE, filename="", overwrite=FALSE, wopt=list())
}

\arguments{
  \item{x}{SpatRaster}
  \item{...}{additional argument of the same type as \code{x} or numeric}
  \item{ties}{character. Indicates how to treat ties. Either "random", "lowest", "highest", "first", or "NA"}
  \item{na.rm}{logical. If \code{TRUE}, \code{NA} values are ignored. If \code{FALSE}, \code{NA} is returned if \code{x} has any \code{NA} values}
  \item{filename}{character. Output filename}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list with named options for writing files as in \code{\link{writeRaster}}}  
}

\value{
SpatRaster
}

\examples{
r <- rast(system.file("ex/logo.tif", package="terra"))   
r <- c(r/2, r, r*2)
m <- modal(r)
}

\keyword{univar}
\keyword{math}
