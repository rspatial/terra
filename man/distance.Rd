\name{distance}

\alias{distance}
\alias{distance,SpatRaster,missing-method}
\alias{distance,SpatRaster,SpatVector-method}

\title{Distance}

\description{
If (\code{y} is missing) this method computes the distance, for all cells that are \code{NA} in SpatRaster \code{x} to the nearest cell that is not \code{NA}.

If (\code{y} is a SpatVector) the distance, for all cells to that SpatVector are computed. For lines and polygons this is done after rasterization (for now).

The distance unit is in meters if the CRS is (\code{+proj=longlat}) and in map units (typically also meters) when it is not.
}

\usage{
\S4method{distance}{SpatRaster,missing}(x, y, filename="", overwrite=FALSE, wopt=list(), ...)
\S4method{distance}{SpatRaster,SpatVector}(x, y, filename="", overwrite=FALSE, wopt=list(), ...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{y}{missing or SpatVector}
  \item{filename}{character. Output filename. Optional}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list. Options for writing files as in \code{\link{writeRaster}}}
  \item{...}{Additional arguments. None implemented}  
}

\value{SpatRaster}


\examples{
r <- rast(ncol=36,nrow=18)
v <- rep(NA, ncell(r))
v[500] <- 1
values(r) <- v
d <- distance(r) 
#plot(d / 100000)

p <- vect(rbind(c(0,0), c(90,30), c(-90,-30)), crs="+proj=longlat +datum=WGS84")
pd <- distance(r, p) 
}

\keyword{spatial}
