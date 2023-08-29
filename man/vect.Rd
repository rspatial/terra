\name{vect}

\docType{methods}

\alias{vect}
\alias{vect,character-method}
\alias{vect,list-method}
\alias{vect,SpatVector-method}
\alias{vect,SpatVectorCollection-method}
\alias{vect,SpatExtent-method}
\alias{vect,data.frame-method}
\alias{vect,matrix-method}
\alias{vect,missing-method}
\alias{vect,Spatial-method}
\alias{vect,sf-method}
\alias{vect,sfc-method}
\alias{vect,XY-method}
\alias{vect,PackedSpatVector-method}


\title{Create SpatVector objects}

\description{
Methods to create a SpatVector from a filename or other R object. 

A filename can be for a shapefile or any spatial file format. 

You can use a data.frame to make a SpatVector of points; or a "geom" matrix to make a SpatVector of any supported geometry (see examples and \code{\link{geom}}). 

You can supply a list of SpatVectors to append them into a single SpatVector. 

SpatVectors can also be created from "Well Known Text", and from spatial vector data objects defined in the \code{sf} or \code{sp} packages.
} 

\usage{
\S4method{vect}{character}(x, layer="", query="", extent=NULL, filter=NULL, 
		crs="", proxy=FALSE, what="")

\S4method{vect}{matrix}(x, type="points", atts=NULL, crs="")

\S4method{vect}{data.frame}(x, geom=c("lon", "lat"), crs="", keepgeom=FALSE)

\S4method{vect}{list}(x, type="points", crs="")

\S4method{vect}{SpatExtent}(x, crs="")

\S4method{vect}{SpatVectorCollection}(x)

\S4method{vect}{sf}(x)
}

\arguments{
\item{x}{character. A filename; or a "Well Known Text" string; SpatExtent, data.frame (to make a SpatVector of points); a "geom" matrix to make a SpatVector of any supported geometry (see examples and \code{\link{geom}}); a spatial vector data object defined in the \code{sf} or \code{sp} packages; or a list with matrices with coordinates}
\item{layer}{character. layer name to select a layer from a file (database) with multiple layers}
\item{query}{character. An query to subset the dataset in the \href{https://gdal.org/user/ogr_sql_dialect.html}{OGR-SQL dialect}}
\item{extent}{Spat* object. The extent of the object is used as a spatial filter to select the geometries to read. Ignored if \code{filter} is not \code{NULL}}
\item{filter}{SpatVector. Used as a spatial filter to select geometries to read (the convex hull is used for lines or points). It is guaranteed that all features that overlap with the extent of filter will be returned. It can happen that additional geometries are returned}
\item{type}{character. Geometry type. Must be "points", "lines", or "polygons"}
\item{atts}{data.frame with the attributes. The number of rows must match the number of geometrical elements}
\item{crs}{character. The coordinate reference system in one of the following formats: WKT/WKT2, <authority>:<code>, or PROJ-string notation (see \code{\link{crs}}). See \code{\link{crs}}}
\item{proxy}{logical. If \code{TRUE} a SpatVectorProxy is returned}
\item{what}{character indicating what to read. Either \code{""} for geometries and attributes, or \code{"geoms"} to only read the geometries, \code{"attributes"} to only read the attributes (that are returned as a data.frame)}
\item{geom}{character. The field name(s) with the geometry data. Either two names for x and y coordinates of points, or a single name for a single column with WKT geometries)}
\item{keepgeom}{logical. If \code{TRUE} the geom variable(s) is (are) also included in the attributes}
}

\seealso{\code{\link{geom}}}


\value{SpatVector}


\examples{
### SpatVector from file
f <- system.file("ex/lux.shp", package="terra")
f
v <- vect(f)
v

## subsetting (large) files
## with attribute query 
v <- vect(f, query="SELECT NAME_1, NAME_2, ID_2 FROM lux WHERE ID_2 < 4")

## with an extent
e <- ext(5.9, 6.3, 49.9, 50) 
v <- vect(f, extent=e)

## with polygons
p <- as.polygons(e)
v <- vect(f, filter=p)


### SpatVector from a geom matrix
x1 <- rbind(c(-180,-20), c(-140,55), c(10, 0), c(-140,-60))
x2 <- rbind(c(-10,0), c(140,60), c(160,0), c(140,-55))
x3 <- rbind(c(-125,0), c(0,60), c(40,5), c(15,-45))
hole <- rbind(c(80,0), c(105,13), c(120,2), c(105,-13))
z <- rbind(cbind(object=1, part=1, x1, hole=0), cbind(object=2, part=1, x3, hole=0),
			cbind(object=3, part=1, x2, hole=0), cbind(object=3, part=1, hole, hole=1))
colnames(z)[3:4] <- c('x', 'y')

p <- vect(z, "polygons")
p

z[z[, "hole"]==1, "object"] <- 4
lns <- vect(z[,1:4], "lines")
plot(p)
lines(lns, col="red", lwd=2)

### from wkt
v <- vect("POLYGON ((0 -5, 10 0, 10 -10, 0 -5))")

wkt <- c("MULTIPOLYGON ( ((40 40, 20 45, 45 30, 40 40)), 
((20 35, 10 30, 10 10, 30 5, 45 20, 20 35),(30 20, 20 15, 20 25, 30 20)))", 
"POLYGON ((0 -5, 10 0, 10 -10, 0 -5))")
w <- vect(wkt)

# combine two SpatVectors 
vw <- rbind(w, v)

# add a data.frame
d <- data.frame(id=1:2, name=c("a", "b"))
values(w) <- d

# add data.frame on creation, here from a geom matrix
g <- geom(w)
d <- data.frame(id=1:2, name=c("a", "b"))
m <- vect(g, "polygons", atts=d, crs="+proj=longlat +datum=WGS84")

### SpatVector from a data.frame
d$wkt <- wkt
x <- vect(d, geom="wkt")

d$wkt <- NULL
d$lon <- c(0,10)
d$lat <- c(0,10)
x <- vect(d, geom=c("lon", "lat"))

# SpatVector to sf
#sf::st_as_sf(x)
}

\keyword{methods}
\keyword{spatial}

