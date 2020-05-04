\name{resample}

\alias{resample}
\alias{resample,SpatRaster,SpatRaster-method}

\title{Transfer values of a SpatRaster to another one with a different geometry}

\description{
resample transfers values between SpatRaster objects that do not align (have a different origin and/or resolution). See \code{\link{project}} to change the coordinate reference system (crs). 

If the origin and crs are the same, you should consider using these other functions instead: \code{\link{aggregate}}, \code{\link{disaggregate}} or \code{\link{crop}}.
}

\usage{
\S4method{resample}{SpatRaster,SpatRaster}(x, y, method="bilinear", 
       filename="", overwrite=FALSE, wopt=list(), ...)
}

\arguments{
  \item{x}{SpatRaster to be resampleed}
  \item{y}{SpatRaster that \code{x} should be resampleed to}
  \item{method}{character. Method used for estimating the new cell values. One of: 
  
	near: nearest neighbour. This method is fast, and it can be the preferred method if the cell values represent classes. It is not a good choice for continuous values.
	bilinear: bilinear interpolation. Default.
	cubic: cubic interpolation.
	cubicspline: cubic spline interpolation.
	lanczos: Lanczos windowed sinc resampling.
  }
  \item{filename}{character. Output filename. Optional}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list. Options for writing files as in \code{\link{writeRaster}}}
  \item{...}{additional arguments. None implemented}  
}


\value{
SpatRaster 
}


\seealso{ \code{\link{aggregate}}, \code{\link{disaggregate}}, \code{\link{crop}}, \code{\link{project}},
}


\examples{

r <- rast(nrow=3, ncol=3, xmin=0, xmax=10, ymin=0, ymax=10)
values(r) <- 1:ncell(r)
s <- rast(nrow=25, ncol=30, xmin=1, xmax=11, ymin=-1, ymax=11)
x <- resample(r, s, method="bilinear")

opar <- par(no.readonly =TRUE)
par(mfrow=c(1,2))
plot(r)
plot(x)
par(opar)
}

\keyword{spatial}

