\name{rasterize}

\docType{methods}

\alias{rasterize}
\alias{rasterize,SpatVector,SpatRaster-method}

\title{rasterize}

\description{
rasterize polygons
}

\usage{
\S4method{rasterize}{SpatVector,SpatRaster}(x, y, background=NA, filename="", ...)

}

\arguments{
  \item{x}{SpatVector object}
  \item{y}{SpatRaster object}
  \item{background}{numeric. Value for cells that are not covered by a polygon}
  \item{filename}{Character, output filename. Optional}
  \item{...}{Additional arguments as for \code{\link{writeRaster}}}  
}



\value{
SpatRaster
}


\examples{
x1 <- rbind(c(-180,-20), c(-140,55), c(10, 0), c(-140,-60))
x2 <- rbind(c(-10,0), c(140,60), c(160,0), c(140,-55))
x3 <- rbind(c(-125,0), c(0,60), c(40,5), c(15,-45))
x4 <- rbind(c(80,0), c(105,13), c(120,2), c(105,-13))

z <- rbind(cbind(object=1, part=1, x1), cbind(object=2, part=1, x2), cbind(object=3, part=1, x3), cbind(object=3, part=2,  x4))
colnames(z)[3:4] <- c('x', 'y')
z <- data.frame(z)
z$hole <- 0
z$hole[z$object==3 & z$part==2] <- 1

p <- SpatPolygon(z)
r <- rast()
a <- rasterize(p, r)


}

\keyword{spatial}
