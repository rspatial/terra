\name{rasterize}

\docType{methods}

\alias{rasterize}
\alias{rasterize,SpatLayer,SpatRaster-method}

\title{rasterize}

\description{
rasterize polygons
}

\usage{
\S4method{rasterize}{SpatLayer,SpatRaster}(x, y, background=NA, filename="", format="", datatype="FLT4S", overwrite=FALSE, ...)

}

\arguments{
  \item{x}{SpatLayer object}
  \item{y}{SpatRaster object}
  \item{background}{numeric. Value for cells that are not covered by a polygon}
  \item{filename}{Character. Output filename. Optional}
  \item{format}{Character. Output file format. See \code{\link{writeRaster}}}
  \item{datatype}{Character. Output file data type. See \code{\link{writeRaster}}}  
  \item{overwrite}{Logical. If \code{TRUE}, \code{filename} is overwritten}
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

