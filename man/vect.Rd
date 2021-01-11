\name{vect}

\docType{methods}

\alias{vect}
\alias{vect,character-method}
\alias{vect,data.frame-method}
\alias{vect,matrix-method}
\alias{vect,missing-method}
\alias{vect,Spatial-method}
\alias{vect,sf-method}
\alias{vect,sfc-method}
\alias{vect,PackedSpatVector-method}


\title{Create SpatVector objects}

\description{
Create a new SpatVector
} 

\usage{
\S4method{vect}{character}(x, ...)

\S4method{vect}{matrix}(x, type="points", atts=NULL, crs="", ...)

\S4method{vect}{data.frame}(x, geom=c("lon", "lat"), crs="", ...)

\S4method{vect}{sf}(x, ...)
}

\arguments{
\item{x}{character (filename or "Well Known Text"); or a data.frame or matrix with geometry data (see \code{\link{geom}}; or missing; or a vector object from \code{sf} or \code{sp}}
\item{type}{character. Geometry type. Must be "points", "lines", or "polygons"}
\item{atts}{data.frame with the attributes. The number of rows must match the number of geometrical elements}
\item{crs}{the coordinate reference system (PROJ4 notation)}
\item{geom}{the field name(s) with the geometry data. Either two names for x and y coordinates of points, or a single name for a single column with WKT geometries)}
\item{...}{additional matrices and/or lists with matrices}  
}

\seealso{\code{\link{geom}}}


\value{SpatVector}


\examples{
### from file
f <- system.file("ex/lux.shp", package="terra")
f
v <- vect(f)
v

### from a geom matrix

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
w <- vect("POLYGON ((0 -5, 10 0, 10 -10, 0 -5))")

wkt <- c("MULTIPOLYGON ( ((40 40, 20 45, 45 30, 40 40)), 
((20 35, 10 30, 10 10, 30 5, 45 20, 20 35),(30 20, 20 15, 20 25, 30 20)))", 
"POLYGON ((0 -5, 10 0, 10 -10, 0 -5))")
w <- vect(wkt)

# add data.frame
g <- geom(w)
d <- data.frame(id=1:2, name=c("a", "b"))
m <- vect(g, atts=d, crs="+proj=longlat +datum=WGS84")

### from a data.frame
d$wkt <- wkt
x <- vect(d, geom="wkt")

d$wkt <- NULL
d$lon <- c(0,10)
d$lat <- c(0,10)
x <- vect(d, geom=c("lon", "lat"))

### SpatVect to sf
#sf::st_as_sf(as.data.frame(w, geom=TRUE), wkt="geometry", crs=crs(w))
}

\keyword{methods}
\keyword{spatial}

