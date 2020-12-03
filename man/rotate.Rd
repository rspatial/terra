\name{rotate}

\docType{methods}

\alias{rotate}
\alias{rotate,SpatRaster-method}
\alias{rotate,SpatVector-method}

\title{Rotate a SpatRaster along longitude}

\description{
Rotate a SpatRaster that has longitude coordinates from 0 to 360, to standard coordinates between -180 and 180 degrees (or vice-versa). Longitude between 0 and 360 is frequently used in global climate models.

Or rotate a SpatVector around the lower-right corner if its extent. 
}

\usage{
\S4method{rotate}{SpatRaster}(x, left=TRUE, filename="", overwrite=FALSE, wopt=list(), ...)

\S4method{rotate}{SpatVector}(x, angle, ...)
}

\arguments{
 \item{x}{SpatRaster}
 \item{left}{logical. If \code{TRUE}, rotate to the left, else to the right}
 \item{filename}{character. Output filename. Optional}
 \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
 \item{wopt}{list. Options for writing files as in \code{\link{writeRaster}}}
 \item{angle}{numeric. Angle of rotation in degrees}
 \item{...}{additional arguments. None implemented}  
}

\value{
SpatRaster
}

\seealso{\code{\link{shift}}}

\examples{
x <- rast(nrow=9, ncol=18, nl=3, xmin=0, xmax=360)
v <- rep(as.vector(t(matrix(1:ncell(x), nrow=9, ncol=18))), 3)
values(x) <- v
z <- rotate(x)
}

\keyword{spatial}
