\name{rasterize}

\docType{methods}

\alias{rasterize}
\alias{rasterize,SpatVector,SpatRaster-method}

\title{Rasterize vector data}

\description{
Transfer vector data to a raster
}

\usage{
\S4method{rasterize}{SpatVector,SpatRaster}(x, y, field=1:nrow(x), background=NA, update=FALSE,
                   touches=is.lines(mask), filename="", overwrite=FALSE, wopt=list(), ...)

}

\arguments{
  \item{x}{SpatVector}
  \item{y}{SpatRaster}
  \item{field}{numeric. The values to be rasterized. Either a field name of \code{x}, a single number, or a vector with the same length as \code{x}}
  \item{background}{numeric. Value to put in the cells that are not covered by any of the features of \code{x}. Default is \code{NA}}
  \item{touches}{logical. If \code{TRUE}, all cells touched by lines or polygons will be updated, not just those on the line render path, or whose center point is within the polygon}
  \item{update}{logical. If \code{TRUE}, the values of the Raster* object are updated for the cells that overlap with the geometries of \code{x}. Default is \code{FALSE}}
  \item{filename}{character. Output filename. Optional}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list. Options for writing files as in \code{\link{writeRaster}}}
  \item{...}{additional arguments. None implemented}  
}

\seealso{
\code{\link{mask}}
}


\value{
SpatRaster
}

\examples{
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
r <- rast(v, ncol=75, nrow=100)
x <- rasterize(v, r)

plot(x)
lines(v)
}


\keyword{spatial}

