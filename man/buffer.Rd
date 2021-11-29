\name{buffer}

\alias{buffer}
\alias{buffer,SpatRaster-method}
\alias{buffer,SpatVector-method}


\title{Create a buffer around vector geometries or raster patches}

\description{
Calculate a buffer around all cells that are not \code{NA} in a SpatRaster, or around the geometries of a SpatVector)
 
Note that the distance unit of the buffer \code{width} parameter is meters if the CRS is (\code{+proj=longlat}), and in map units (typically also meters) if not. 
}

\usage{
\S4method{buffer}{SpatRaster}(x, width, filename="", ...)

\S4method{buffer}{SpatVector}(x, width, quadsegs=10)
}

\arguments{
\item{x}{SpatRaster or SpatVector}
\item{width}{numeric. Unit is meter if \code{x} has a longitude/latitude CRS, or mapunits in other cases. Should be > 0 for SpatRaster}
  \item{filename}{character. Output filename}
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
  \item{quadsegs}{positive integer. Number of line segments to use to draw a quart circle}  
}


\value{SpatRaster}

\seealso{
\code{\link{distance}}
}


\examples{
r <- rast(ncols=36, nrows=18)
v <- rep(NA, ncell(r))
v[500] <- 1
values(r) <- v
b <- buffer(r, width=5000000) 
plot(b)

v <- vect(rbind(c(10,10), c(0,60)), crs="+proj=merc")
b <- buffer(v, 20)
plot(b)
points(v)

crs(v) <- "+proj=longlat" 
b <- buffer(v, 1500000)
plot(b)
points(v)
}

\keyword{spatial}
