\name{subset}

\alias{subset}
\alias{subset,SpatRaster-method}
\alias{[[,SpatRaster,numeric,missing-method}
\alias{[[,SpatRaster,logical,missing-method}
\alias{[[,SpatRaster,character,missing-method}
\alias{[,SpatRaster,character,missing-method}
\alias{[,SpatRasterStack,numeric,missing-method}
\alias{[,SpatRasterStack,numeric,numeric-method}
\alias{[,SpatRasterStack,logical,missing-method}
\alias{[,SpatRasterStack,character,missing-method}
\alias{$,SpatRaster-method}
\alias{$,SpatRasterStack-method}

\title{Subset of a SpatRaster}

\description{
Select a subset of layers from a SpatRaster.  
}

\usage{
\S4method{subset}{SpatRaster}(x, subset, filename="", overwrite=FALSE, wopt=list(), ...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{subset}{integer or character. Should indicate the layers (represented as integer or by their names)}
  \item{filename}{character. Output filename. Optional}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list. Options for writing files as in \code{\link{writeRaster}}}
  \item{...}{additional arguments. None implemented}  
}

\value{
SpatRaster
}

\examples{
s <- rast(system.file("ex/logo.tif", package="terra"))   
subset(s, 2:3)
subset(s, c(3,2,3,1))
#equivalent to 
s[[ c(3,2,3,1) ]]

s[c("red", "green")]
s$red
}

\keyword{ spatial }
