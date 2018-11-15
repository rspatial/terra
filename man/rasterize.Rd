\name{rasterize}

\docType{methods}

\alias{rasterize}
\alias{rasterize,SpatLayer,SpatRaster-method}

\title{rasterize}

\description{
rasterize polygons
}

\usage{
\S4method{rasterize}{SpatLayer,SpatRaster}(x, y, background=NA, filename="", overwrite=FALSE, wopt=list(), ...)

}

\arguments{
  \item{x}{SpatLayer}
  \item{y}{SpatRaster}
  \item{background}{numeric. Value for cells that are not covered by a polygon}
  \item{filename}{character. Output filename. Optional}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list. Options for writing files as in \code{\link{writeRaster}}}
  \item{...}{Additional arguments. None implemented}  
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

