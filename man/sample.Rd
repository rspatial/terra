\name{spatSample}

\docType{methods}

\alias{spatSample}
\alias{spatSample,SpatRaster-method}

\title{Take a regular sample}

\description{
Take a regular sample of a SpatRaster. Either get cell values, or (when \code{as.raster=TRUE}) get a new SpatRaster with the same extent, but fewer cells. Note that, in order to assure regularity when requesting a regular sample, the number of values returned may not be exactly the same as the \code{size} requested.
}

\usage{
\S4method{spatSample}{SpatRaster}(x, size, method="regular", replace=FALSE, as.raster=FALSE, cells=FALSE, ...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{size}{numeric. The sample size}
  \item{method}{character. Should be "regular" or "random"}
  \item{replace}{logical. If \code{TRUE}, sampling is with replacement (if \code{method="random"}}
  \item{as.raster}{logical. If \code{TRUE}, a SpatRaster is returned}
  \item{cells}{logical. If \code{TRUE}, cellnumbers are returned instead of values}
  \item{...}{additional arguments. None implemented}
}


\details{
In some cases you may want to know the location of the sampled cells. In that situation you can take a sample of the cell numbers and use extract. See examples.
}

\value{
numeric or SpatRaster
}

\examples{
f <- system.file("ex/test.tif", package="terra")
r <- rast(f)
s <- spatSample(r, 10, as.raster=TRUE)
spatSample(r, 10)
spatSample(r, 10, "random")

## if you require cell numbers and/or coordinates
size <- 6
# random cells 
cells <- spatSample(r, 6, "random", cells=TRUE)
v <- r[cells]
xy <- xyFromCell(r, cells)
cbind(xy, v)

# regular 
cells <- spatSample(r, 6, "regular", cells=TRUE)
v <- r[cells]
xy <- xyFromCell(r, cells)
cbind(xy, v)
}

\keyword{spatial}
