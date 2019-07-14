\name{collapse}

\alias{collapse}
\alias{collapse,SpatRaster-method}


\title{Collapse cell values from a multi-layer SpatRaster}

\description{
Use a single layer SpatRaster object to select cell values from different layers in a multi-layer SpatRaster and "collapse" it into a single layer. The values of the SpatRaster to select layers (\code{y}) should be between \code{1} and \code{nlyr(x)} (values outside this range are ignored); they are also truncated to integers. 

See \code{\link{extract}} for extraction of values by cell, point, or otherwise.
}

\usage{
\S4method{collapse}{SpatRaster}(x, y, filename="", overwrite=FALSE, wopt=list(), ...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{y}{SpatRaster}
  \item{filename}{character. Output filename. Optional}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list. Options for writing files as in \code{\link{writeRaster}}}
  \item{...}{additional arguments. None implemented}  
}

\value{
SpatRaster
}

\seealso{\code{\link{tapp}}, \code{\link{extract}}}


\examples{
r <- rast(ncol=10, nrow=10)
values(r) <- 1
s <- c(r, r+2, r+5)
set.seed(1)
values(r) <- round((runif(ncell(r)))*3)
x <- collapse(s, r)
}

\keyword{methods}
\keyword{spatial}

