\name{extract}

\docType{methods}

\alias{extract}
\alias{extract,SpatRaster,SpatVector-method}
\alias{extract,SpatRaster,SpatExtent-method}
\alias{extract,SpatRaster,matrix-method}
\alias{extract,SpatRaster,data.frame-method}
\alias{extract,SpatRaster,numeric-method}
\alias{extract,SpatVector,SpatVector-method}


\alias{[,SpatVector,SpatVector,missing-method}
\alias{[,SpatVector,SpatExtent,missing-method}
\alias{[,SpatRaster,SpatVector,missing-method}
\alias{[,SpatRaster,missing,missing-method}
\alias{[,SpatRaster,numeric,missing-method}
\alias{[,SpatRaster,missing,numeric-method}
\alias{[,SpatRaster,numeric,numeric-method}
\alias{[,SpatRaster,SpatRaster,missing-method}
\alias{[,SpatRaster,SpatExtent,missing-method}


\title{Extract values from a SpatRaster}

\description{
Extract values from a SpatRaster for a set of locations. The locations can be a SpatVector (points, lines, polygons), a matrix with (x, y) or (longitude, latitude -- in that order!) coordinates, or a vector with cell numbers.  

When argument \code{y} is a \code{SpatVector}, and \code{list=FALSE}, the first column has the ID (record number) of the \code{SpatVector} used.
} 

\usage{
\S4method{extract}{SpatRaster,SpatVector}(x, y, fun=NULL, method="simple", list=FALSE, factors=TRUE, 
        cells=FALSE, xy=FALSE, weights=FALSE, exact=FALSE,
        touches=is.lines(y), layer=NULL, ...)

\S4method{extract}{SpatRaster,SpatExtent}(x, y, factors=TRUE, cells=FALSE, xy=FALSE)

\S4method{extract}{SpatRaster,matrix}(x, y, ...)

\S4method{extract}{SpatRaster,numeric}(x, y, ...)
}

\arguments{
\item{x}{SpatRaster}
\item{y}{SpatVector (for points, lines, polygons), or for points, 2-column matrix or data.frame (x, y) or (lon, lat), or a vector with cell numbers}
\item{fun}{function to summarize the data by geometry. If \code{weights=TRUE} or \code{exact=TRUE} only code{mean}, \code{sum}, \code{min} and \code{max} are accepted)}
\item{...}{additional arguments to \code{fun} if \code{y} is a SpatVector. For example \code{na.rm=TRUE}. Or arguments passed to the \code{SpatRaster,SpatVector} method if \code{y} is a matrix (such as the \code{method} and \code{cells} arguments)}
\item{method}{character. method for extracting values with points ("simple" or "bilinear"). With "simple" values for the cell a point falls in are returned. With "bilinear" the returned values are interpolated from the values of the four nearest raster cells}
\item{list}{logical. If \code{FALSE} the output is simplified to a \code{matrix} (if \code{fun=NULL})}
\item{factors}{logical. If \code{TRUE} the categories are returned as factors instead of their numerical representation. The value returned becomes a data.frame if it otherwise would have been a matrix, even if there are no factors}
\item{cells}{logical. If \code{TRUE} the cell numbers are also returned, unless \code{fun} is not \code{NULL}. Also see \code{\link{cells}}}
\item{xy}{logical. If \code{TRUE} the coordinates of the cells are also returned, unless \code{fun} is not \code{NULL}. Also see \code{\link{xyFromCell}}}
\item{weights}{logical. If \code{TRUE} and \code{y} has polygons, the approximate fraction of each cell that is covered is returned as well, for example to compute a weighted mean}
\item{exact}{logical. If \code{TRUE} and \code{y} has polygons, the exact fraction of each cell that is covered is returned as well, for example to compute a weighted mean}
\item{touches}{logical. If \code{TRUE}, values for all cells touched by lines or polygons are extracted, not just those on the line render path, or whose center point is within the polygon. Not relevant for points; and always considered \code{TRUE} when \code{weights=TRUE} or \code{exact=TRUE}}
\item{layer}{character or numeric to select the layer to exctract from for each geometry. If \code{layer} is a character it can be a name in \code{y} or a vector of layer names. If it is numeric, it must be integer values between \code{1} and \code{nlyr(x)}}
}

\value{matrix, list, or data.frame}

\seealso{\link{values}}

\examples{
r <- rast(ncols=5, nrows=5, xmin=0, xmax=5, ymin=0, ymax=5)
values(r) <- 1:25
xy <- rbind(c(0.5,0.5), c(2.5,2.5))
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
z <- rast(v, resolution=.1, names="test")
values(z) <- 1:ncell(z)

rf <- system.file("ex/elev.tif", package="terra")
x <- rast(rf)
#extract(x, v, mean, na.rm=TRUE)

#e <- extract(z, v)
#e
#tapply(e[,2], e[,1], mean, na.rm=TRUE)

#ee <- extract(z, v, list=TRUE)
#rapply(ee, mean)

#x <- c(z, z*2, z/3)
#names(x) <- letters[1:3]


#e <- extract(x, v)
#de <- data.frame(e)
#aggregate(de[,2:4], de[,1,drop=FALSE], mean)

#ee <- extract(x, v, list=TRUE)
#matrix(rapply(ee, mean), ncol=nlyr(x), byrow=TRUE)
}

\keyword{methods}
\keyword{spatial}

