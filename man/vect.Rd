\name{vect}

\docType{methods}

\alias{vect}
\alias{vect,character-method}
\alias{vect,data.frame-method}
\alias{vect,matrix-method}
\alias{vect,missing-method}
\alias{vect,Spatial-method}
\alias{vect,PackedSpatVector-method}

\title{Create SpatVector objects}

\description{
Create a new SpatVector
} 

\usage{
\S4method{vect}{data.frame}(x, type="points", atts=NULL, crs=NA, ...)
}

\arguments{
\item{x}{character (filename), data.frame, matrix, or missing. If \code{x} is a filename, all other arguments are ignored}
\item{type}{character. Geometry type. Must be "points", "lines", or "polygons"}
\item{atts}{data.frame with the attributes. The number of rows must match the number of geometrical elements}
\item{crs}{the coordinate reference system (PROJ4 notation)}
\item{...}{additional matrices and/or lists with matrices}  
}


\value{SpatVector object}


\examples{

f <- system.file("exdata/lux.shp", package="terra")
f
v <- vect(f)
v

x1 <- rbind(c(-180,-20), c(-140,55), c(10, 0), c(-140,-60))
x2 <- rbind(c(-10,0), c(140,60), c(160,0), c(140,-55))
x3 <- rbind(c(-125,0), c(0,60), c(40,5), c(15,-45))
hole <- rbind(c(80,0), c(105,13), c(120,2), c(105,-13))

z <- rbind(cbind(object=1, part=1, x1, hole=0), cbind(object=2, part=1, x3, hole=0),
			cbind(object=3, part=1, x2, hole=0), cbind(object=3, part=1, hole, hole=1))
colnames(z)[3:4] <- c('x', 'y')
z <- data.frame(z)

p <- vect(z, "polygons")
p

z[z$hole==1, "object"] <- 4
lns <- vect(z[,1:4], "lines")
plot(p)
lines(lns, col='red')

}

\keyword{methods}
\keyword{spatial}

