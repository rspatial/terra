\name{extract}

\docType{methods}

\alias{extract}
\alias{extract,SpatRaster,SpatVector-method}
\alias{extract,SpatRaster,sf-method}
\alias{extract,SpatRaster,SpatExtent-method}
\alias{extract,SpatRaster,matrix-method}
\alias{extract,SpatRaster,data.frame-method}
\alias{extract,SpatRaster,numeric-method}
\alias{extract,SpatVector,SpatVector-method}
\alias{extract,SpatVector,matrix-method}
\alias{extract,SpatVector,data.frame-method}
\alias{extract,SpatRasterCollection,ANY-method}
\alias{extract,SpatRasterDataset,ANY-method}


\title{Extract values from a SpatRaster}

\description{
Extract values from a SpatRaster for a set of locations. The locations can be a SpatVector (points, lines, polygons), a data.frame or matrix with (x, y) or (longitude, latitude -- in that order!) coordinates, or a vector with cell numbers.  

When argument \code{y} is a \code{SpatVector} the first column has the ID (record number) of the \code{SpatVector} used (unless you set \code{ID=FALSE}).

Alternatively, you can use \code{\link{zonal}} after using \code{\link{rasterize}} with a \code{SpatVector} (this may be more efficient in some cases).
} 

\usage{
\S4method{extract}{SpatRaster,SpatVector}(x, y, fun=NULL, method="simple", cells=FALSE, xy=FALSE,
    ID=TRUE, weights=FALSE, exact=FALSE, touches=is.lines(y), small=TRUE,
	layer=NULL, bind=FALSE, raw=FALSE, search_radius=0, ...)

\S4method{extract}{SpatRaster,SpatExtent}(x, y, cells=FALSE, xy=FALSE)

\S4method{extract}{SpatRaster,matrix}(x, y, cells=FALSE, method="simple")

\S4method{extract}{SpatRaster,numeric}(x, y, xy=FALSE, raw=FALSE)

\S4method{extract}{SpatVector,SpatVector}(x, y, count=FALSE)
}

\arguments{
\item{x}{SpatRaster or SpatVector of polygons}
\item{y}{SpatVector (points, lines, or polygons). Alternatively, for points, a 2-column matrix or data.frame (x, y) or (lon, lat). Or a vector with cell numbers}
\item{fun}{function to summarize the extracted data by line or polygon geometry. You can use \code{fun=table} to tabulate raster values for each line or polygon geometry. If \code{weights=TRUE} or \code{exact=TRUE} only \code{mean}, \code{sum}, \code{min}, \code{max} and \code{table} are accepted). Ignored if \code{y} has point geometry}
\item{method}{character. method for extracting values with points ("simple" or "bilinear"). With "simple" values for the cell a point falls in are returned. With "bilinear" the returned values are interpolated from the values of the four nearest raster cells}
\item{cells}{logical. If \code{TRUE} the cell numbers are also returned, unless \code{fun} is not \code{NULL}. Also see \code{\link{cells}}}
\item{xy}{logical. If \code{TRUE} the coordinates of the cells are also returned, unless \code{fun} is not \code{NULL}. See \code{\link{xyFromCell}}}
\item{ID}{logical. Should an ID column be added? If so, the first column returned has the IDs (record numbers) of \code{y}}
\item{weights}{logical. If \code{TRUE} and \code{y} has polygons, the approximate fraction of each cell that is covered is returned as well, for example to compute a weighted mean}
\item{exact}{logical. If \code{TRUE} and \code{y} has polygons, the exact fraction of each cell that is covered is returned as well, for example to compute a weighted mean}
\item{touches}{logical. If \code{TRUE}, values for all cells touched by lines or polygons are extracted, not just those on the line render path, or whose center point is within the polygon. Not relevant for points; and always considered \code{TRUE} when \code{weights=TRUE} or \code{exact=TRUE}}
\item{small}{logical. If \code{TRUE}, values for all cells in touched polygons are extracted if none of the cells center points is within the polygon; even if \code{touches=FALSE}}
\item{layer}{character or numeric to select the layer to extract from for each geometry. If \code{layer} is a character it can be a name in \code{y} or a vector of layer names. If it is numeric, it must be integer values between \code{1} and \code{nlyr(x)}}
\item{bind}{logical. If \code{TRUE}, a SpatVector is returned consisting of the input SpatVector \code{y} and the \code{cbind}-ed extracted values}
\item{raw}{logical. If \code{TRUE}, a matrix is returned with the "raw" numeric cell values. If \code{FALSE}, a data.frame is returned and the cell values are transformed to factor, logical, or integer values, where appropriate}
\item{search_radius}{positive number. A search-radius that is used when \code{y} has point geometry. If this value is larger than zero, it is the maximum distance used to find the a cell with a value that is nearest to the cell that the point falls in if that cell that has a missing (\code{NA}) value. The value of this nearest cell, the distance to the original cell, and the new cell number are returned. The radius should be expressed in m if the data have lon/lat coordinates or in the distance unit of the crs in other cases (typically also m). For lon/lat data, the mean latitude of the points is used to compute the distances, so this may be imprecise for data with a large latitudinal range}
\item{...}{additional arguments to \code{fun} if \code{y} is a SpatVector. For example \code{na.rm=TRUE}. Or arguments passed to the \code{SpatRaster,SpatVector} method if \code{y} is a matrix (such as the \code{method} and \code{cells} arguments)}
\item{count}{logical. If \code{TRUE} and \code{x} has polygons geometry and \code{y} has points geometry, the number of points in polygons is returned}
}

\value{data.frame, matrix or SpatVector}

\seealso{\code{\link{values}, \link{zonal}, \link{extractAlong}, \link{extractRange}, \link{rapp}}}

\examples{
r <- rast(ncols=5, nrows=5, xmin=0, xmax=5, ymin=0, ymax=5)
values(r) <- 1:25
xy <- cbind(lon=c(0.5,2.5), lat=c(0.5,2.5))
p <- vect(xy, crs="+proj=longlat +datum=WGS84")

extract(r, xy)
extract(r, p)

r[1,]
r[5]
r[,5]

r[c(0:2, 99:101)]

f <- system.file("ex/meuse.tif", package="terra")
r <- rast(f)

xy <- cbind(179000, 330000)
xy <- rbind(xy-100, xy, xy+1000)
extract(r, xy)

p <- vect(xy)
g <- geom(p)
g

extract(r, p)

x <- r + 10
extract(x, p)

i <- cellFromXY(r, xy)
x[i]
r[i]

y <- c(x,x*2,x*3)
y[i]

## extract with a polygon
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
v <- v[1:2,]

rf <- system.file("ex/elev.tif", package="terra")
x <- rast(rf)
extract(x, v, mean, na.rm=TRUE)

z <- rast(v, resolution=.1, names="test")
values(z) <- 1:ncell(z)
e <- extract(z, v, ID=TRUE)
e
tapply(e[,2], e[,1], mean, na.rm=TRUE)

x <- c(z, z*2, z/3)
names(x) <- letters[1:3]

e <- extract(x, v, ID=TRUE)
de <- data.frame(e)
aggregate(de[,2:4], de[,1,drop=FALSE], mean)
}

\keyword{methods}
\keyword{spatial}

