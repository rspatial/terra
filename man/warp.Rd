\name{warp}

\alias{warp}
\alias{warp,SpatRaster,SpatRaster-method}

\title{Transfer values of a SpatRaster to another one with a different geometry}

\description{
Warp transfers values between SpatRaster objects that do not align (have a different origin and/or coordinate reference system (crs)). This is also known as "resample" and/or "reproject". If the origin and crs are the same, you should consider using these other functions instead: \code{\link{aggregate}}, \code{\link{disaggregate}}, \code{\link{crop}}.
}

\usage{
\S4method{warp}{SpatRaster,SpatRaster}(x, y, method="bilinear", 
       filename="", overwrite=FALSE, wopt=list(), ...)
}

\arguments{
  \item{x}{SpatRaster to be warped}
  \item{y}{SpatRaster that \code{x} should be warped to}
  \item{method}{character. Method used to compute values. \code{"bilinear"} for bilinear interpolation, or \code{"ngb"} for "nearest neighbor"}
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
x <- warp(r, s, method='bilinear')

par(mfrow=c(1,2))
plot(r)
plot(x)

# warp to different CRS
a <- rast(ncol=40, nrow=40, xmin=-110, xmax=-90, ymin=40, ymax=60, 
          crs="+proj=longlat +datum=WGS84")
values(a) = 1:ncell(a)
b <- rast(ncol=94, nrow=124, xmin=-944881, xmax=935118, ymin=4664377, ymax=7144377, 
     crs="+proj=lcc +lat_1=48 +lat_2=33 +lon_0=-100 +datum=WGS84")
w <- warp(a, b)
}

\keyword{spatial}

