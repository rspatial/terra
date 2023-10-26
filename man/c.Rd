\name{c}

\docType{methods}

\alias{c}
\alias{c,SpatRaster-method}
\alias{c,SpatRasterDataset-method}
\alias{c,SpatRasterCollection-method}
\alias{c,SpatVector-method}
\alias{c,SpatVectorCollection-method}


\title{Combine SpatRaster or SpatVector objects}

\description{
With \code{c} you can:

-- Combine \code{SpatRaster} objects. They must have the same extent and resolution. However, if \code{x} is empty (has no cell values), its geometry is ignored with a warning. Two empty SpatRasters with the same geometry can also be combined (to get a summed number of layers). Also see \code{\link{add<-}}

-- Add a \code{SpatRaster} to a \code{SpatRasterDataset} or \code{SpatRasterCollection}

-- Add \code{SpatVector} objects to a new or existing \code{SpatVectorCollection}

To append SpatVectors, use \code{rbind}. 
}

\seealso{\code{\link{add<-}}}


\usage{
\S4method{c}{SpatRaster}(x, ..., warn=TRUE)

\S4method{c}{SpatRasterDataset}(x, ...)

\S4method{c}{SpatRasterCollection}(x, ...)

\S4method{c}{SpatVector}(x, ...)

\S4method{c}{SpatVectorCollection}(x, ...)
}

\arguments{
  \item{x}{SpatRaster, SpatVector, SpatRasterDataset or SpatVectorCollection}
  \item{warn}{logical. If \code{TRUE}, a warning is emitted if \code{x} is an empty SpatRaster}
  \item{...}{as for \code{x} (you can only combine raster with raster data and vector with vector data)}
}


\value{
Same class as \code{x}
}


\examples{
r <- rast(nrows=5, ncols=9)
values(r) <- 1:ncell(r)
x <- c(r, r*2, r*3)
}

\keyword{spatial}

