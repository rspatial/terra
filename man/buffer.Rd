\name{buffer}

\alias{buffer}
\alias{buffer,SpatRaster-method}
\alias{buffer,SpatVector-method}


\title{Create a buffer around vector objects or raster patches}

\description{
Calculate a buffer around all cells that are not \code{NA} in a SpatRaster, or around the objcts in a SpatVector (currently only implemented for points)
 
Note that the distance unit of the buffer \code{width} parameter is meters if the CRS is (\code{+proj=longlat}), and in map units (typically also meters) if not. 
}

\usage{
\S4method{buffer}{SpatRaster}(x, width, filename="", overwrite=FALSE, wopt=list(), ...)
\S4method{buffer}{SpatVector}(x, width, quadsegs=10, capstyle="round", ...)
}

\arguments{
\item{x}{SpatRaster}
\item{width}{numeric. Unit is meter if \code{x} has a longitude/latitude CRS, or mapunits in other cases. Should be > 0 for SpatRaster}
  \item{filename}{character. Output filename. Optional}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list. Options for writing files as in \code{\link{writeRaster}}}
  \item{...}{additional arguments. None implemented}  
  \item{quadsegs}{postive integer. Number of line segments to use to draw a quart circle}  
  \item{capstyle}{character. Style of cap to use at the ends of the geometry. Allowed values: "round", "flat", "square" (ignored for now) }  
}


\value{SpatRaster}

\seealso{
\code{\link{distance}}
}


\examples{
r <- rast(ncol=36,nrow=18)
v <- rep(NA, ncell(r))
v[500] <- 1
values(r) <- v
b <- buffer(r, width=5000000) 
plot(b)

v <- vect(rbind(c(10,10), c(0,60)))
b <- buffer(v, 20)
plot(b)
points(v)

crs(v) <- "+proj=longlat +datum=WGS84"
b <- buffer(v, 1500000)
plot(b)
points(v)
}

\keyword{spatial}
