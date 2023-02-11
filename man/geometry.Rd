\name{geom}

\docType{methods}

\alias{geom}
\alias{geom,SpatVector-method}

\title{Get the geometry (coordinates) of a SpatVector}

\description{
Get the geometry of a SpatVector. If \code{wkt=FALSE}, this is a five-column matrix or data.frame: the vector object ID, the IDs for the parts of each object (e.g. five polygons that together are one spatial object), the x (longitude) and y (latitude) coordinates, and a flag indicating whether the part is a "hole" (only relevant for polygons). 

If \code{wkt=TRUE}, the "well-known text" representation is returned as a character vector.
}

\usage{
\S4method{geom}{SpatVector}(x, wkt=FALSE, hex=FALSE, df=FALSE, list=FALSE, xnm="x", ynm="y")
}

\arguments{
  \item{x}{SpatVector}
  \item{wkt}{logical. If \code{TRUE} the WKT geometry is returned (unless \code{hex} is also \code{TRUE})}  
  \item{hex}{logical. If \code{TRUE} the hexadecimal geometry is returned}  
  \item{df}{logical. If \code{TRUE} a \code{data.frame} is returned instead of a matrix (only if \code{wkt=FALSE}, \code{hex=FALSE}, and \code{list=FALSE})} 
  \item{list}{logical. If \code{TRUE} a nested \code{list} is returned with data.frames of coordinates}
  \item{xnm}{character. If \code{list=TRUE} the "x" column name for the coordinates data.frame }  
  \item{ynm}{character. If \code{list=TRUE} the "y" column name for the coordinates data.frame}  
}


\value{
matrix, vector, data.frame, or list
}

\seealso{
\code{\link{crds}}, \code{\link{xyFromCell}}
}


\examples{
x1 <- rbind(c(-175,-20), c(-140,55), c(10, 0), c(-140,-60))
x2 <- rbind(c(-125,0), c(0,60), c(40,5), c(15,-45))
x3 <- rbind(c(-10,0), c(140,60), c(160,0), c(140,-55))
x4 <- rbind(c(80,0), c(105,13), c(120,2), c(105,-13))
z <- rbind(cbind(object=1, part=1, x1), cbind(object=2, part=1, x2), 
           cbind(object=3, part=1, x3), cbind(object=3, part=2,  x4))
colnames(z)[3:4] <- c('x', 'y')
z <- cbind(z, hole=0)
z[(z[, "object"]==3 & z[,"part"]==2), "hole"] <- 1

p <- vect(z, "polygons")
geom(p)

f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
g <- geom(v)
head(g)

w <- geom(v, wkt=TRUE)
substr(w, 1, 60)
}

\keyword{spatial}
