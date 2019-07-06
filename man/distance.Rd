\name{distance}

\alias{distance}
\alias{distance,SpatRaster,missing-method}
\alias{distance,SpatRaster,SpatVector-method}
\alias{distance,SpatVector,missing-method}
\alias{distance,SpatVector,SpatVector-method}

\title{Geographic distance}

\description{

\bold{If \code{x} is a SpatRaster:}

If \code{y} is \code{missing} this method computes the distance, for all cells that are \code{NA} in SpatRaster \code{x} to the nearest cell that is not \code{NA}.

If \code{y} is a SpatVector, the distance to that SpatVector is computed for all cells. For lines and polygons this is done after rasterization (for now).

\bold{If \code{x} is a SpatVector:}

If \code{y} is \code{missing}, a distance matrix between all object in \code{x} is computed. An distance matrix object of class "dist" is returned. 

If \code{y} is a SpatVector the geographic distance between all elements two sets of objects is computed (and a matrix is returned). If both sets have the same number of points, and \code{pairwise=TRUE}, the distance between each pair of objects is computed, and a vector is returned.
}

\usage{
\S4method{distance}{SpatRaster,missing}(x, y, filename="", overwrite=FALSE, wopt=list(), ...)

\S4method{distance}{SpatRaster,SpatVector}(x, y, filename="", overwrite=FALSE, wopt=list(), ...)

\S4method{distance}{SpatVector,missing}(x, y, ...)

\S4method{distance}{SpatVector,SpatVector}(x, y, pairwise=FALSE, ...)
}

\arguments{
  \item{x}{SpatRaster or SpatVector}
  \item{y}{missing or SpatVector}
  \item{filename}{character. Output filename. Optional}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list. Options for writing files as in \code{\link{writeRaster}}}
  \item{...}{Additional arguments. None implemented}  
  \item{pairwise}{logical. If \code{TRUE} and if x and y have the same size (number of rows), the pairwise distances are returned instead of the distances between all elements} 
}

\value{
SpatRaster or numeric or matrix or distance matrix (object of class "dist")

The unit is in meters if the CRS is \code{longlat} and in map units (typically also meters) when it is not.
}


\examples{
r <- rast(ncol=36,nrow=18)
v <- rep(NA, ncell(r))
v[500] <- 1
values(r) <- v
d <- distance(r) 
#plot(d / 100000)

p1 <- vect(rbind(c(0,0), c(90,30), c(-90,-30)), crs="+proj=longlat +datum=WGS84")
dp <- distance(r, p1) 


d <- distance(p1)
d
as.matrix(d)

p2 <- vect(rbind(c(30,-30), c(25,40), c(-9,-3)), crs="+proj=longlat +datum=WGS84")
dd <- distance(p1, p2)
dd
pd <- distance(p1, p2, pairwise=TRUE)
pd
pd == diag(dd)
}

\keyword{spatial}

