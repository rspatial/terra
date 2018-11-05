\name{rasterize}

\docType{methods}

\alias{rasterize}
\alias{rasterize,SpatLayer,SpatRaster-method}

\title{rasterize}

\description{
rasterize polygons
}

\usage{
\S4method{rasterize}{SpatLayer,SpatRaster}(x, y, background=NA, filename="", ...)

}

\arguments{
  \item{x}{SpatLayer object}
  \item{y}{SpatRaster object}
  \item{background}{numeric. Value for cells that are not covered by a polygon}
  \item{filename}{Character, output filename. Optional}
  \item{...}{Additional arguments as for \code{\link{writeRaster}}}  
}


\value{
SpatRaster
}

\examples{
f <- system.file("exdata/lux.shp", package="terra")
v <- vect(f)
r <- rast(v, ncol=75, nrow=100)
x <- rasterize(v, r)

image(x)
lines(v)
}


\keyword{spatial}

