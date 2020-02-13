\name{spatSample}

\docType{methods}

\alias{spatSample}
\alias{spatSample,SpatRaster-method}

\title{Take a regular sample}

\description{
Take a regular sample of a SpatRaster. Either get cell values, or a new SpatRaster with the same extent, but fewer cells. 
}

\usage{
\S4method{spatSample}{SpatRaster}(x, size, method="regular", replace=FALSE, as.raster=FALSE, ...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{size}{numeric. The sample size}
  \item{method}{character. Should be "regular" or "random"}
  \item{replace}{logical. If \code{TRUE}, sampling is with replacement (if \code{method="random"}}
  \item{as.raster}{logical. If \code{TRUE} and \code{method="regular"}, a SpatRaster is returned}
  \item{...}{additional arguments. None implemented}
}



\value{
numeric or SpatRaster
}


\examples{
f <- system.file("exdata/test.tif", package="terra")
r <- rast(f)
s <- spatSample(r, 10, as.raster=TRUE)
spatSample(r, 10)
spatSample(r, 10, "random")
}

\keyword{spatial}
