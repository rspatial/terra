\name{gridDistance}

\alias{gridDistance}
\alias{gridDistance,SpatRaster-method}


\title{Distance on a grid}

\description{
The function calculates the distance to cells of a SpatRaster when the path has to go through the centers of the eight neighboring raster cells.
 
The distance is in meters if the coordinate reference system (CRS) of the SpatRaster is longitude/latitude (\code{+proj=longlat}) and in the units of the CRS (typically meters) in other cases. 
 
Distances are computed by summing local distances between cells, which are connected with their neighbors in 8 directions.

The shortest distance to the cells with the \code{target} value is computed for all cells that are not \code{NA}. Cells that are \code{NA} cannot be traversed and are ignored, unless the target itself is \code{NA}, in which case the distance to the nearest cell that is not \code{NA} is computed for all cells that are \code{NA}. 
}

\usage{
\S4method{gridDistance}{SpatRaster}(x, target=0, scale=1000, maxiter=50, filename="", ...) 
}

\arguments{
\item{x}{SpatRaster}
\item{target}{numeric. value of the target cells (where to compute distance to)}
\item{scale}{numeric. Scale factor for longitude/latitude data (1 = m, 1000 = km)}
\item{maxiter}{numeric. The maximum number of iterations. Increase this number if you get the warning that \code{costDistance did not converge}}
\item{filename}{character. output filename (optional)}
\item{...}{additional arguments as for \code{\link{writeRaster}}}  
}


\seealso{See \code{\link[terra]{distance}} for "as the crow flies" distance} 


\note{
see terra:::.oldGridDistance for the previous version of this method (will be removed in a future version)
}

\value{SpatRaster}



\examples{
# global lon/lat raster
r <- rast(ncol=10,nrow=10, vals=1)
r[48] <- 0
r[66:68] <- NA
d <- gridDistance(r) 
plot(d)


# planar
crs(r) <- "+proj=utm +zone=15 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
d <- gridDistance(r) 
plot(d)

# distance to cells that are not NA 
rr <- classify(r, cbind(1, NA))
dd <- gridDistance(rr, NA) 


}

\keyword{spatial}
