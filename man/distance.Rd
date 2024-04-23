\name{distance}

\alias{distance}
\alias{distance,SpatRaster,missing-method}
\alias{distance,SpatRaster,SpatVector-method}
\alias{distance,SpatRaster,sf-method}
\alias{distance,SpatVector,ANY-method}
\alias{distance,SpatVector,SpatVector-method}
\alias{distance,matrix,matrix-method}
\alias{distance,matrix,missing-method}
\alias{distance,data.frame,data.frame-method}
\alias{distance,data.frame,missing-method}

\title{Geographic distance}

\description{

If \code{x} is a \bold{SpatRaster}:

If \code{y} is \code{missing} this method computes the distance, for all cells that are \code{NA} in SpatRaster \code{x} to the nearest cell that is not \code{NA} (or other values, see arguments "target" and "exclude"). 

If \code{y} is a numeric value, the cells with that value are ignored. That is, distance to or from these cells is not computed (only if \code{grid=FALSE}).

If \code{y} is a SpatVector, the distance to that SpatVector is computed for all cells. For lines and polygons this is done after rasterization; and only the overlapping areas of the vector and raster are considered (for now).

The distance is always expressed in meter if the coordinate reference system is longitude/latitude, and in map units otherwise. Map units are typically meter, but inspect \code{crs(x)} if in doubt.

Results are more precise, sometimes much more precise, when using longitude/latitude rather than a planar coordinate reference system, as these distort distance. 

If \code{x} is a \bold{SpatVector}:

If \code{y} is \code{missing}, a distance matrix between all objects in \code{x} is computed. A distance matrix object of class "dist" is returned. 

If \code{y} is a SpatVector the geographic distance between all objects is computed (and a matrix is returned). If both sets have the same number of points, and \code{pairwise=TRUE}, the distance between each pair of objects is computed, and a vector is returned.

The distance is always expressed in meter, except when the coordinate reference system is longitude/latitude AND one of the SpatVector(s) consists of lines or polygons. In that case the distance is in degrees, and thus not very useful (this will be fixed soon). Otherwise, results are more precise, sometimes much more precise, when using longitude/latitude rather than a planar coordinate reference system, as these distort distance. 


If \code{x} is a \bold{matrix}:

\code{x} should consist of two columns, the first with "x" (or longitude) and the second with "y" coordinates (or latitude). If \code{y} is a also a matrix, the distance between each points in \code{x} and all points in \code{y} is computed, unless \code{pairwise=TRUE}

If \code{y} is missing, the distance between each points in \code{x} with all other points in \code{x} is computed, unless \code{sequential=TRUE}
}



\usage{
\S4method{distance}{SpatRaster,missing}(x, y, target=NA, exclude=NULL, unit="m", haversine=TRUE, filename="", ...)

\S4method{distance}{SpatRaster,SpatVector}(x, y, unit="m", rasterize=FALSE, haversine=TRUE, filename="", ...)

\S4method{distance}{SpatVector,ANY}(x, y, sequential=FALSE, pairs=FALSE, symmetrical=TRUE, unit="m")

\S4method{distance}{SpatVector,SpatVector}(x, y, pairwise=FALSE, unit="m")

\S4method{distance}{matrix,matrix}(x, y, lonlat, pairwise=FALSE, unit="m")

\S4method{distance}{matrix,missing}(x, y, lonlat, sequential=FALSE, pairs=FALSE, symmetrical=TRUE, unit="m")
}

\arguments{
  \item{x}{SpatRaster, SpatVector, or two-column matrix with coordinates (x,y) or (lon,lat)}
  \item{y}{missing, numeric, SpatVector, or two-column matrix}
  \item{target}{numeric. The value of the cells for which distances to cells that are not \code{NA} should be computed}
  \item{exclude}{numeric. The value of the cells that should not be considered for computing distances}
  \item{unit}{character. Can be either "m" or "km"}
  \item{haversine}{logical. Use the haversine formula for lon/lat data use the haversine formula? If \code{FALSE}, the more precise but slower method of Karney (2003) is used}  
  \item{rasterize}{logical. If \code{TRUE} distance is computed from the cells covered by the geometries after rasterization. This can be much faster in some cases}
  \item{filename}{character. Output filename}
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
  \item{sequential}{logical. If \code{TRUE}, the distance between sequential geometries is returned}
  \item{pairwise}{logical. If \code{TRUE} and if x and y have the same size (number of rows), the pairwise distances are returned instead of the distances between all elements} 
  \item{lonlat}{logical. If \code{TRUE} the coordinates are interpreted as angular (longitude/latitude). If \code{FALSE} they are interpreted as planar}
  \item{pairs}{logical. If \code{TRUE} a "from", "to", "distance" matrix is returned} 
  \item{symmetrical}{logical. If \code{TRUE} and \code{pairs=TRUE}, the distance between a pair is only included once. The distance between geometry 1 and 3 is included, but the (same) distance between 3 and 1 is not}
}

\value{
SpatRaster or numeric or matrix or distance matrix (object of class "dist")
}

\note{
A distance matrix can be coerced into a regular matrix with \code{as.matrix}
}

\references{
Karney, C.F.F., 2013. Algorithms for geodesics, J. Geodesy 87: 43-55. doi:10.1007/s00190-012-0578-z.
}

\examples{
#lonlat
r <- rast(ncols=36, nrows=18, crs="+proj=longlat +datum=WGS84")
r[500] <- 1
d <- distance(r) 
plot(d / 100000)

#planar
rr <- rast(ncols=36, nrows=18, crs="+proj=utm +zone=1 +datum=WGS84")
rr[500] <- 1
d <- distance(rr) 

rr[3:10, 3:10] <- 99
e <- distance(rr, exclude=99) 

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


# polygons, lines
crs <- "+proj=utm +zone=1"
p1 <- vect("POLYGON ((0 0, 8 0, 8 9, 0 9, 0 0))", crs=crs)
p2 <- vect("POLYGON ((5 6, 15 6, 15 15, 5 15, 5 6))", crs=crs)
p3 <- vect("POLYGON ((2 12, 3 12, 3 13, 2 13, 2 12))", crs=crs)
p <- rbind(p1, p2, p3)
L1 <- vect("LINESTRING(1 11, 4 6, 10 6)", crs=crs)
L2 <- vect("LINESTRING(8 14, 12 10)", crs=crs)
L3 <- vect("LINESTRING(1 8, 12 14)", crs=crs)
lns <- rbind(L1, L2, L3)
pts <- vect(cbind(c(7,10,10), c(3,5,6)), crs=crs)

distance(p1,p3)
distance(p)
distance(p,pts)
distance(p,lns)
distance(pts,lns)
}

\keyword{spatial}

