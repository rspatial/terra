\name{crop}

\docType{methods}

\alias{crop}
\alias{crop,SpatRaster-method}

\title{Crop}

\description{
crop returns a new SpatRaster that is the geographic subset of an object as specified by a SpatExtent object.
}

\usage{
\S4method{crop}{SpatRaster}(x, y, snap='near', filename="", overwrite=FALSE, ...)

}

\arguments{
  \item{x}{SpatRaster object}
  \item{y}{SpatExtent object}
  \item{snap}{Character. One of 'near', 'in', or 'out'}  
  \item{filename}{Character, output filename. Optional}
  \item{overwrite}{Logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{...}{Additional arguments as for \code{\link{writeRaster}}}  
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
