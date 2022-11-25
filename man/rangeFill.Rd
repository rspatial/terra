\name{rangeFill}

\alias{rangeFill}
\alias{rangeFill,SpatRaster-method}

\title{Fill layers with a range}

\description{
Fill layers with cell-varying ranges defined by a start and end SpatRaster. 
The range must start at 1 and end at a user-defined maximum. 
}

\usage{
\S4method{rangeFill}{SpatRaster}(x, y, limit, filename="", ...) 
}

\arguments{
\item{x}{SpatRaster with at single layer. Cell values are indicate the start of the range (1 based)}

\item{x}{SpatRaster with at single layer. Cell values are indicate the end of the range}

\item{limit}{numeric > 1. The range size}

\item{filename}{character. Output filename}

\item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}

\item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}

}

\value{SpatRaster}


\seealso{ \code{\link{rapp}} } 


\examples{

library(terra)
x <- y <- rast(ncol=2, nrow=2)
values(x) <- c(NA, 1:3)
values(y) <- c(NA, 4:6)
r <- rangeFill(x, y, 8)

}


\keyword{spatial}

