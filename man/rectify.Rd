\name{rectify}

\alias{rectify}
\alias{rectify,SpatRaster-method}


\title{rectify a SpatRaster}


\description{
Rectify a rotated SpatRaster into a non-rotated object
}

\usage{
\S4method{rectify}{SpatRaster}(x, method="bilinear", ext=NULL, filename="", overwrite=FALSE, wopt=list(), ...)
}

\arguments{
  \item{x}{SpatRaster to be rectified}
  \item{method}{character. Method used to for resampling. See \code{\link{resample}}}
  \item{ext}{SpatExtent to crop \code{x} to a smaller area of interest}
  \item{filename}{character. Output filename. Optional}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list. Options for writing files as in \code{\link{writeRaster}}}
  \item{...}{additional arguments. None implemented}    
}


\value{
SpatRaster
}


\keyword{spatial}

