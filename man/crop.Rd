\name{crop}

\docType{methods}

\alias{crop}
\alias{crop,SpatRaster-method}

\title{Crop}

\description{
Cut out a part of a SpatRaster with a SpatExtent, or another object from which an extent can be obtained. You can only extract rectangular areas, but see \code{\link{mask}} for setting cell values within SpatRaster to \code{NA}.
}

\usage{
\S4method{crop}{SpatRaster}(x, y, snap="near", filename="", overwrite=FALSE, ...)

}

\arguments{
  \item{x}{SpatRaster object}
  \item{y}{SpatExtent object}
  \item{snap}{Character. One of "near", "in", or "out"}  
  \item{filename}{Character. Output filename. Optional}
  \item{overwrite}{Logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{...}{Additional arguments. None implemented}  
}



\value{
SpatRaster
}


\examples{
r <- rast(nrow=45, ncol=90)
values(r) <- 1:ncell(r)
e <- ext(-160, 10, 30, 60)
rc <- crop(r, e)
}

\keyword{spatial}
