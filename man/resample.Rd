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
    
	\code{bilinear}: bilinear interpolation (3x3 cell window). This is used by default if the first layer of \code{x} is not categorical

	\code{average}:  This can be a good choice with continuous variables if the output cells overlap with multiple input cells. 

	\code{near}: nearest neighbor. This is used by default if the first layer of \code{x} is categorical. This method is not a good choice for continuous values.

	\code{mode}: The modal value. This can be a good choice for categrical rasters, if the output cells overlap with multiple input cells. 
  
	\code{cubic}: cubic interpolation (5x5 cell window).

	\code{cubicspline}: cubic B-spline interpolation. (5x5 cell window).

	\code{lanczos}: Lanczos windowed sinc resampling. (7x7 cell window).

	\code{sum}: the weighted sum of all non-NA contributing grid cells.

	\code{min, q1, median, q3, max}: the minimum, first quartile, median, third quartile, or maximum value. 
	
	\code{rms}: the root-mean-square value of all non-NA contributing grid cells.
  }
  
  
  \item{threads}{logical. If \code{TRUE} multiple threads are used (faster for large files)}
  
  \item{by_util}{logical. If \code{TRUE} the GDAL warp utility is used}
 
  \item{filename}{character. Output filename}

  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
}


\value{SpatRaster }


\seealso{ \code{\link{aggregate}}, \code{\link{disagg}}, \code{\link{crop}}, \code{\link{project}}}


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

