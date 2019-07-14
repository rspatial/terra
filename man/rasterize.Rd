\name{rasterize}

\docType{methods}

\alias{rasterize}
\alias{rasterize,SpatVector,SpatRaster-method}

\title{Rasterize vector data}

\description{
Transfer vector data to a raster
}

\usage{
\S4method{rasterize}{SpatVector,SpatRaster}(x, y, background=NA, filename="", overwrite=FALSE, wopt=list(), ...)

}

\arguments{
  \item{x}{SpatVector}
  \item{y}{SpatRaster}
  \item{background}{numeric. Value for cells that are not covered by a polygon}
  \item{filename}{character. Output filename. Optional}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list. Options for writing files as in \code{\link{writeRaster}}}
  \item{...}{additional arguments. None implemented}  
}


\value{
SpatRaster
}

\examples{
f <- system.file("exdata/lux.shp", package="terra")
v <- vect(f)
r <- rast(v, ncol=75, nrow=100)
x <- rasterize(v, r)

plot(x)
lines(v)
}


\keyword{spatial}

