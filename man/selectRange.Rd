\name{selectRange}

\alias{selectRange}
\alias{selectRange,SpatRaster-method}


\title{Select the values of a range of layers, as specified by cell values in another SpatRaster}

\description{
Use a single layer SpatRaster object to select cell values from different layers in a multi-layer SpatRaster. The values of the SpatRaster to select layers (\code{y}) should be between \code{1} and \code{nlyr(x)} (values outside this range are ignored); they are also truncated to integers. 

See \code{\link{rapp}} for applying af function to a range of variable size.

See \code{\link{extract}} for extraction of values by cell, point, or otherwise.
}

\usage{
\S4method{selectRange}{SpatRaster}(x, y, z=1, repint=0, filename="", overwrite=FALSE, wopt=list(), ...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{y}{SpatRaster. Cell values must be positive integers. They indicate the first layer to select for each cell}
  \item{z}{positive integer. The number of layers to select}
  \item{repint}{integer > 1 and < nlyr(x) allowing for repeated selection at a fixed interval. For example, if \code{x} has 36 layers, and the value of a cell in \code{y}=2 and \code{repint} = 12, the values for layers 2, 14 and 26 are returned}
  \item{filename}{character. Output filename. Optional}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list. Options for writing files as in \code{\link{writeRaster}}}
  \item{...}{additional arguments. None implemented}  
}

\value{
SpatRaster
}

\seealso{\code{\link{rapp}}, \code{\link{tapp}}, \code{\link{extract}}}


\examples{
r <- rast(ncol=10, nrow=10)
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

