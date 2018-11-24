\name{warp}

\alias{warp}
\alias{warp,SpatRaster,SpatRaster-method}

\title{warp a SpatRaster object}

\description{
Warp transfers values between SpatRaster objects that do not align (have a different origin). This is also known as "resample".

Adjusting to different coordinate reference system (projection) is not implemented yet.

Before using warp, you may want to consider using these other functions instead: \code{\link[raster]{aggregate}}, \code{\link[raster]{disaggregate}}, \code{\link[raster]{crop}}, and others.
}


\usage{
\S4method{warp}{SpatRaster,SpatRaster}(x, y, method="bilinear", filename="", overwrite=FALSE, wopt=list(), ...)
}

\arguments{
  \item{x}{SpatRaster to be warped}
  \item{y}{SpatRaster that \code{x} should be warped to}
  \item{method}{character. Method used to compute values. \code{"bilinear"} for bilinear interpolation, or \code{"ngb"} for using the nearest neighbor}
  \item{filename}{character. Output filename. Optional}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list. Options for writing files as in \code{\link{writeRaster}}}
  \item{...}{Additional arguments. None implemented}  
}


\value{
SpatRaster 
}


\seealso{ \code{\link[raster]{aggregate}}, \code{\link[raster]{disaggregate}}, \code{\link[raster]{crop}} }


\examples{

r <- rast(nrow=3, ncol=3, xmin=0, xmax=10, ymin=0, ymax=10)
values(r) <- 1:ncell(r)
s <- rast(nrow=25, ncol=30, xmin=1, xmax=11, ymin=-1, ymax=11)
x <- warp(r, s, method='bilinear')
#par(mfrow=c(1,2))
#plot(r)
#plot(s)
}

\keyword{spatial}

