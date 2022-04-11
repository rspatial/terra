\name{costDistance}

\alias{costDistance}
\alias{costDistance,SpatRaster-method}


\title{Cost distance}

\description{

Compute the cost-distance to cells that have a value of zero on a cost surface, and assuming that the path has to go through the centers of one of the neighboring raster cells. 

Distances are computed by summing local distances between cells, which are connected with their neighbors in 8 directions.
 
The cost surface must express the cost per unit distance (speed) of travel. 

If the coordinate reference system (CRS) of the SpatRaster is longitude/latitude (\code{+proj=longlat}) the distance is in meters divided by variable m. The default of m is 1000, expressing distance in km. Otherwise the distance is in the units of the CRS (typically meters).

}

\usage{
\S4method{costDistance}{SpatRaster}(x, scale=1000, filename="", ...) 
}

\arguments{
\item{x}{SpatRaster}
\item{scale}{scale factor for longlat data (1 = m, 1000 = km)}
\item{filename}{character. output filename (optional)}
\item{...}{additional arguments as for \code{\link{writeRaster}}}  
}


\seealso{\code{\link{gridDistance}, \link[terra]{distance}} } 


\value{SpatRaster}


\examples{
r <- rast(ncols=5, nrows=5, crs="+proj=utm +zone=1 +datum=WGS84", xmin=0, xmax=5, ymin=0, ymax=5, vals=1)
r[13] <- 0
d <- costDistance(r)
plot(d)
text(d, digits=1)
}

\keyword{spatial}
