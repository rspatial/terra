\name{cover}

\docType{methods}

\alias{cover}
\alias{cover,SpatRaster,SpatRaster-method}


\title{Cover (replace) \code{NA} values with values of another raster }

\description{
Replace \code{NA} values (or another value) in SpatRaster (\code{x}) with the values of SpatRaster (\code{y})
}

\usage{
\S4method{cover}{SpatRaster,SpatRaster}(x, y, value=NA, filename="", wopt=list(), ...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{y}{SpatRaster}
  \item{value}{numeric. The cell values in \code{x} to be replaced by the values in \code{y}}

  \item{filename}{character. Output filename. Optional}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list. Options for writing files as in \code{\link{writeRaster}}}
  \item{...}{Additional arguments. None implemented}  
}


\value{
SpatRaster
}

\examples{
r1 <- r2 <- rast(ncols=36, nrows=18)
values(r1) <- 1:ncell(r1)
values(r2) <- runif(ncell(r2))
r2 <- reclassify(r2, cbind(-Inf, 0.5, NA))
r3 <- cover(r2, r1)
}


\keyword{methods}
\keyword{spatial}
