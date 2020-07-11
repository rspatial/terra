\name{crop}

\docType{methods}

\alias{crop}
\alias{crop,SpatRaster,ANY-method}
\alias{crop,SpatVector,SpatVector-method}
\alias{crop,SpatVector,SpatExtent-method}

\title{Cut out a geographic subset}

\description{
Cut out a part of a SpatRaster with a SpatExtent, or another object from which an extent can be obtained. You can only extract rectangular areas, but see \code{\link{mask}} for setting cell values within SpatRaster to \code{NA}.
}

\usage{
\S4method{crop}{SpatRaster,ANY}(x, y, snap="near", filename="", overwrite=FALSE, wopt=list(), ...)

\S4method{crop}{SpatVector,SpatVector}(x, y, ...)

}

\arguments{
  \item{x}{SpatRaster or SpatVector}
  \item{y}{SpatExtent or other object that has a SpatExtent; or SpatVector if \code{x} is a SpatVector}
  \item{snap}{character. One of "near", "in", or "out"}  
  \item{filename}{character. Output filename. Optional}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list. Options for writing files as in \code{\link{writeRaster}}}
  \item{...}{additional arguments. None implemented}  
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

