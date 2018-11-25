\name{crop}

\docType{methods}

\alias{crop}
\alias{crop,SpatRaster-method}

\title{Crop}

\description{
Cut out a part of a SpatRaster with a SpatExtent, or another object from which an extent can be obtained. You can only extract rectangular areas, but see \code{\link{mask}} for setting cell values within SpatRaster to \code{NA}.
}

\usage{
\S4method{crop}{SpatRaster}(x, y, snap="near", filename="", overwrite=FALSE, wopt=list(), ...)

}

\arguments{
  \item{x}{SpatRaster}
  \item{y}{SpatExtent}
  \item{snap}{Character. One of "near", "in", or "out"}  
  \item{filename}{character. Output filename. Optional}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list. Options for writing files as in \code{\link{writeRaster}}}
  \item{...}{Additional arguments. None implemented}  
}



\value{
SpatRaster
}


\examples{
r <- rast(xmin=0, xmax=10, ymin=0, ymax=10, nrows=25, ncols=25)
values(r) <- 1:ncell(r)
e <- ext(-5, 5, -5, 5)
rc <- crop(r, e)
}

\keyword{spatial}
