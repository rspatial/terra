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
p1 <- rbind(c(-180,-20), c(-140,55), c(10, 0), c(-140,-60), c(-180,-20))
p2 <- rbind(c(-10,0), c(140,60), c(160,0), c(140,-55), c(-10,0))
p3 <- rbind(c(-125,0), c(0,60), c(40,5), c(15,-45), c(-125,0))
x <- spPolygons(list(p1, p2), p3)
r <- georst()

a <- rasterize(x, r)
}

\keyword{spatial}
