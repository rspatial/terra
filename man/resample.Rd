\name{resample}

\alias{resample}
\alias{resample,SpatRaster,SpatRaster-method}

\title{Transfer values of a SpatRaster to another one with a different geometry}

\description{
resample transfers values between SpatRaster objects that do not align (have a different origin and/or resolution). See \code{\link{project}} to change the coordinate reference system (crs). 

If the origin and extent of the input and output are the same, you should consider using these other functions instead: \code{\link{aggregate}}, \code{\link{disagg}}, \code{\link{extend}} or \code{\link{crop}}.
}

\usage{
\S4method{resample}{SpatRaster,SpatRaster}(x, y, method, threads=FALSE, filename="", ...)
}

\arguments{
  \item{x}{SpatRaster to be resampled}
  \item{y}{SpatRaster with the geometry that \code{x} should be resampled to}
  \item{method}{character. Method used for estimating the new cell values. One of: 
  
	\code{near}: nearest neighbor. This method is fast, and it can be the preferred method if the cell values represent classes. It is not a good choice for continuous values. This is used by default if the first layer of \code{x} is categorical.

	\code{bilinear}: bilinear interpolation. This is the default if the first layer of \code{x} is numeric (not categorical).

	\code{cubic}: cubic interpolation.

	\code{cubicspline}: cubic spline interpolation.

	\code{lanczos}: Lanczos windowed sinc resampling.

	\code{sum}: the weighted sum of all non-NA contributing grid cells.

	\code{min, q1, med, q3, max, average, mode, rms}: the minimum, first quartile, median, third quartile, maximum, mean, mode, or root-mean-square value of all non-NA contributing grid cells.

  }
  \item{threads}{logical. If \code{TRUE} multiple threads are used (faster for large files)}
  
  \item{filename}{character. Output filename}

  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
}


\value{
SpatRaster 
}


\seealso{ \code{\link{aggregate}}, \code{\link{disagg}}, \code{\link{crop}}, \code{\link{project}}
}


\examples{

r <- rast(nrows=3, ncols=3, xmin=0, xmax=10, ymin=0, ymax=10)
values(r) <- 1:ncell(r)
s <- rast(nrows=25, ncols=30, xmin=1, xmax=11, ymin=-1, ymax=11)
x <- resample(r, s, method="bilinear")

opar <- par(no.readonly =TRUE)
par(mfrow=c(1,2))
plot(r)
plot(x)
par(opar)
}

\keyword{spatial}

