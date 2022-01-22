\name{gridDistance}

\alias{gridDistance}
\alias{gridDistance,SpatRaster-method}


\title{Distance on a grid}

\description{
The function calculates the distance to cells of a SpatRaster when the path has to go through the centers of the eight neighboring raster cells.
 
The distance is in meters if the coordinate reference system (CRS) of the SpatRaster is longitude/latitude (\code{+proj=longlat}) and in the units of the CRS (typically meters) in other cases. 
 
Distances are computed by summing local distances between cells, which are connected with their neighbors in 8 directions.
}

\usage{
\S4method{gridDistance}{SpatRaster}(x, origin, omit=NULL, chunk=FALSE, filename="", overwrite=FALSE, ...) 
}

\arguments{
\item{x}{SpatRaster}
\item{origin}{value(s) of the cells from which the distance is calculated. If \code{origin=NULL} all cells that are not \code{NA} are origins}
\item{omit}{value(s) of the cells that cannot be traversed (optional)}
\item{chunk}{logical. If \code{TRUE} the data are processed in chunks. This may be necessary for large datasets, but could lead to errors in the case of complex areas that need to be omitted if they are spread over different chunks (meandering rivers, for instance).}
\item{filename}{character. output filename (optional)}
\item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
\item{...}{additional arguments as for \code{\link{writeRaster}}}  
}


\seealso{See \code{\link[terra]{distance}} for "as the crow flies" distance} 


\details{
If the SpatRaster to be processed is big, it will be processed in chunks. This may lead to errors in the case of complex objects that cannot be traversed are spread over different chunks (meandering rivers, for instance). 


Only the first layer of SpatRaster \code{x} is processed. See \code{\link{sapp}} or \code{\link{lapply}} to iterate over all layers
}

\value{SpatRaster}


\author{Jacob van Etten and Robert J. Hijmans }


\examples{

#world lon/lat raster
r <- rast(ncol=10,nrow=10, vals=1)
r[48] <- 2
r[66:68] <- 3
d <- gridDistance(r,origin=2,omit=3) 
plot(d)

#UTM small area
crs(r) <- "+proj=utm +zone=15 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
d <- gridDistance(r,origin=2,omit=3) 
plot(d)

}

\keyword{spatial}
