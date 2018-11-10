\name{subset}

\alias{subset}
\alias{subset,SpatRaster-method}
\alias{[[,SpatRaster,numeric,missing-method}
\alias{$,SpatRaster-method}

\title{Subset a SpatRaster object}

\description{
Subset a set of layers from a SpatRaster object.  
}

\usage{
\S4method{subset}{SpatRaster}(x, subset, filename="", format="", datatype="FLT4S", overwrite=FALSE, ...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{subset}{integer or character. Should indicate the layers (represented as integer or by their names)}
  \item{filename}{Character. Output filename. Optional}
  \item{format}{Character. Output file format. See \code{\link{writeRaster}}}
  \item{datatype}{Character. Output file data type. See \code{\link{writeRaster}}} 
  \item{overwrite}{Logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{...}{Additional arguments. None implemented}
}

\value{
Raster* object
}

\examples{
s <- rast(system.file("exdata/logo.grd", package="terra"))   
subset(s, 2:3)
subset(s, c(3,2,3,1))
#equivalent to 
s[[ c(3,2,3,1) ]]
}

\keyword{ spatial }
