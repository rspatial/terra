\name{subset}

\alias{subset}
\alias{subset,SpatRaster-method}
\alias{[[,SpatRaster,numeric,missing-method}
\alias{[[,SpatRaster,logical,missing-method}
\alias{[[,SpatRaster,character,missing-method}
\alias{[,SpatRaster,character,missing-method}
\alias{[,SpatRasterDataset,numeric,missing-method}
\alias{[,SpatRasterDataset,numeric,numeric-method}
\alias{[,SpatRasterDataset,logical,missing-method}
\alias{[,SpatRasterDataset,character,missing-method}
\alias{[[,SpatRasterDataset,ANY,ANY-method}

\alias{[,SpatRasterCollection,numeric,missing-method}

\alias{$,SpatRaster-method}
\alias{$,SpatRasterDataset-method}

\title{Subset of a SpatRaster}

\description{
Select a subset of layers from a SpatRaster.  
}

\usage{
\S4method{subset}{SpatRaster}(x, subset, filename="", overwrite=FALSE, ...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{subset}{integer or character. Should indicate the layers (represented as integer or by their names)}
  \item{filename}{character. Output filename}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
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

s[[c("red", "green")]]
s$red

# expresion based (partial) matching of names with single brackets
s["re"]
s["^re"]

# not with double brackets
# s[["re"]]

}

\keyword{ spatial }
