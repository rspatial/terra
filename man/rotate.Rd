\name{rotate}

\docType{methods}

\alias{rotate}
\alias{rotate,SpatRaster-method}

\title{rotate}

\description{
Rotate a SpatRaster that has x coordinates (longitude) from 0 to 360, to standard coordinates between -180 and 180 degrees (or vice/versa). Longitude between 0 and 360 is frequently used in global climate models.
}

\usage{
\S4method{rotate}{SpatRaster}(x, left=TRUE, filename="", overwrite=FALSE, wopt=list(), ...)
}

\arguments{
 \item{x}{Raster* object}
 \item{left}{logical. If \code{TRUE}, rotate to the left, else to the right}
 \item{filename}{character. Output filename. Optional}
 \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
 \item{wopt}{list. Options for writing files as in \code{\link{writeRaster}}}
 \item{...}{additional arguments. None implemented}  
}

\value{
SpatRaster
}

\examples{
x <- rast(nrow=9, ncol=18, xmin=0, xmax=360)
v <- as.vector(t(matrix(1:ncell(x), nrow=9, ncol=18)))
values(x) <- v
z <- rotate(x)
}

\keyword{spatial}
