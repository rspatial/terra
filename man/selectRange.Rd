\name{selectRange}

\alias{selectRange}
\alias{selectRange,SpatRaster-method}


\title{Select the values of a range of layers, as specified by cell values in another SpatRaster}

\description{
Use a single layer SpatRaster to select cell values from different layers in a multi-layer SpatRaster. The values of the SpatRaster to select layers (\code{y}) should be whole numbers between \code{1} and \code{nlyr(x)} (values outside this range are ignored). 

See \code{\link{rapp}} for applying a function to a range of variable size.

See \code{\link{extract}} for extraction of values by cell, point, or otherwise.
}

\usage{
\S4method{selectRange}{SpatRaster}(x, y, z=1, repint=0, filename="", ...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{y}{SpatRaster. Cell values must be positive integers. They indicate the first layer to select for each cell}
  \item{z}{positive integer. The number of layers to select}
  \item{repint}{integer > 1 and < nlyr(x) allowing for repeated selection at a fixed interval. For example, if \code{x} has 36 layers, and the value of a cell in \code{y}=2 and \code{repint} = 12, the values for layers 2, 14 and 26 are returned}
  \item{filename}{character. Output filename}
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
}

\value{
SpatRaster
}

\seealso{\code{\link{rapp}}, \code{\link{tapp}}, \code{\link{extract}}}


\examples{
r <- rast(ncols=10, nrows=10)
values(r) <- 1
s <- c(r, r+2, r+5)
s <- c(s, s)
set.seed(1)
values(r) <- sample(3, ncell(r), replace=TRUE)
x <- selectRange(s, r)

x <- selectRange(s, r, 3)
}

\keyword{methods}
\keyword{spatial}

