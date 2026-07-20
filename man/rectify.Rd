\name{rectify}

\alias{rectify}
\alias{rectify,SpatRaster-method}


\title{Rectify a SpatRaster}


\description{
Rectify a rotated SpatRaster into a non-rotated object
}

\usage{
\S4method{rectify}{SpatRaster}(x, method="bilinear", aoi=NULL, snap=TRUE,
         threads=FALSE, filename="", ...)
}

\arguments{
  \item{x}{SpatRaster to be rectified}
  \item{method}{character. Method used to for resampling. See \code{\link{resample}}}
  \item{aoi}{SpatExtent or SpatRaster to crop \code{x} to a smaller area of interest; Using a SpatRaster allowing to set the exact output extent and output resolution}
  \item{snap}{logical. If \code{TRUE}, the origin and resolution of the output are the same as would the case when \code{aoi = NULL}. Only relevant if \code{aoi} is a \code{SpatExtent}}
  \item{threads}{logical or positive integer. If \code{TRUE}, multiple threads are used (faster for large files), no more than the "threads" setting of \code{\link{terraOptions}}. A number sets the thread count for this call directly}
  \item{filename}{character. Output filename}
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
}


\seealso{\code{\link{is.rotated}}}

\value{
SpatRaster
}


\keyword{spatial}

