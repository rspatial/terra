\name{rectify}

\alias{rectify}


\title{rectify a SpatRaster}


\description{
Rectify a rotated SpatRaster into a non-rotated object
}

\usage{
rectify(x, method="bilinear",  filename="", overwrite=FALSE, wopt=list())
}

\arguments{
  \item{x}{Raster* object to be rectified}
  \item{method}{character. Method used to for resampling. See \code{\link{resample}}}
  \item{filename}{character. Output filename. Optional}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list. Options for writing files as in \code{\link{writeRaster}}}
}


\value{
SpatRaster
}


\keyword{spatial}

