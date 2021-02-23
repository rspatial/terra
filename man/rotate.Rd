\name{rotate}

\docType{methods}

\alias{rotate}
\alias{rotate,SpatRaster-method}

\title{Rotate a SpatRaster along longitude}

\description{
Rotate a SpatRaster that has longitude coordinates from 0 to 360, to standard coordinates between -180 and 180 degrees (or vice-versa). Longitude between 0 and 360 is frequently used in global climate models.
}

\usage{
\S4method{rotate}{SpatRaster}(x, left=TRUE, filename="", ...)
}

\arguments{
 \item{x}{SpatRaster or SpatVector}
 \item{left}{logical. If \code{TRUE}, rotate to the left, else to the right}
 \item{filename}{character. Output filename}
 \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
}

\value{
SpatRaster
}

\seealso{\code{\link{shift}} and \code{\link{spin}}}

\examples{
x <- rast(nrow=9, ncol=18, nl=3, xmin=0, xmax=360)
v <- rep(as.vector(t(matrix(1:ncell(x), nrow=9, ncol=18))), 3)
values(x) <- v
z <- rotate(x)
}

\keyword{spatial}
