\name{coords}

\docType{methods}

\alias{coords}
\alias{coords,SpatVector-method}
\alias{coords,SpatRaster-method}

\title{Get the coordinates of SpatVector geometries or SpatRaster cells}

\description{
Get the coordinates of a SpatVector or SpatRaster cells. A matrix or data.frame of the x (longitude) and y (latitude) coordinates is returned. 
}

\usage{
\S4method{coords}{SpatVector}(x, df=FALSE)

\S4method{coords}{SpatRaster}(x, df=FALSE)
}

\arguments{
  \item{x}{SpatVector}
  \item{df}{logical. If \code{TRUE} a \code{data.frame} is returned in stead of a matrix}  
}


\value{
matrix or data.frame
}

\seealso{
\code{\link{geom}} returns the complete structure of SpatVector geometries
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
coords(p)

f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
g <- coords(v)
head(g)

}

\keyword{spatial}
