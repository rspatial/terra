\name{gridDistance}

\alias{gridDistance}
\alias{gridDistance,SpatRaster-method}


\title{Distance on a grid}

\description{
The function calculates the distance to cells on a raster when the path has to go through the centers of neighboring cells. Distances are computed by summing local distances between cells, which are connected with their neighbours in 8 directions ('queen' case). 

In future versions, the distance is expressed in meters if the coordinate reference system (CRS) is longitude/latitude (\code{+proj=longlat}) and in the units of the CRS (typically meters) in other cases. 
 
}

\usage{
\S4method{gridDistance}{SpatRaster}(x, filename="", overwrite=FALSE, wopt=list(), ...) 
}

\arguments{
  \item{x}{SpatRaster}
  \item{filename}{character. Output filename. Optional}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list. Options for writing files as in \code{\link{writeRaster}}}
  \item{...}{Additional arguments. None implemented}  

}



\value{SpatRaster}



\examples{
r <- rast(ncol=10, nrow=10)
v <- rep(NA, 100)
v[23] <- 1
values(r) <- v
g <- gridDistance(r)
}

\keyword{spatial}
