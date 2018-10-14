\name{subset}

\alias{subset}
\alias{subset,SpatRaster-method}
\alias{[[,SpatRaster,numeric,missing-method}

\title{Subset a SpatRaster object}

\description{
Subset a set of layers from a SpatRaster object.  
}

\usage{
\S4method{subset}{SpatRaster}(x, subset, filename='', overwrite=FALSE, ...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{subset}{integer or character. Should indicate the layers (represented as integer or by their names)}
  \item{filename}{character. Output filename (optional)}
  \item{overwrite}{Logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{...}{additional arguments as for \code{\link{writeRaster}}}
}

\value{
Raster* object
}

\examples{
s <- rast(system.file("external/rlogo.grd", package="raster"))   
subset(s, 2:3)
subset(s, c(3,2,3,1))
#equivalent to 
s[[ c(3,2,3,1) ]]
}

\keyword{ spatial }
