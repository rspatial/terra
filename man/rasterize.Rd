\name{rasterize}

\docType{methods}

\alias{rasterize}
\alias{rasterize,SpatVector,SpatRaster-method}

\title{Rasterize vector data}

\description{
Transfer vector data to a raster
}

\usage{
\S4method{rasterize}{SpatVector,SpatRaster}(x, y, field, fun, background=NA, update=FALSE,
                   touches=is.lines(x), filename="", overwrite=FALSE, wopt=list(), ...)
}

\arguments{
  \item{x}{SpatVector}
  \item{y}{SpatRaster}
  \item{field}{character, numeric, or missing. If \code{field} is a character, it should a variable name in \code{x} or \code{fun}, or vector of values if \code{fun} can handle these. If \code{field} is numeric it should be a single number, or a vector with the same length as \code{x}. If it is missing, \code{1:nrow(x)} is used}
  \item{fun}{function}{summarizing function that returns a single number. For when there are multiple points in one cell. For example \code{mean}, \code{length} (to get a count), \code{min} or \code{max}. Only used if \code{x} consists of points}
  \item{background}{numeric. Value to put in the cells that are not covered by any of the features of \code{x}. Default is \code{NA}}
  \item{touches}{logical. If \code{TRUE}, all cells touched by lines or polygons will be updated, not just those on the line render path, or whose center point is within the polygon}
  \item{update}{logical. If \code{TRUE}, the values of the SpatRaster are updated for the cells that overlap with the geometries of \code{x}. Default is \code{FALSE}}
  \item{filename}{character. Output filename. Optional}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list. Options for writing files as in \code{\link{writeRaster}}}
  \item{...}{additional arguments passed on to \code{fun}}  
}

\seealso{
\code{\link{mask}}
}


\value{
SpatRaster
}

\examples{
#f <- system.file("ex/lux.shp", package="terra")
#v <- vect(f)
#r <- rast(v, ncol=75, nrow=100)
#x <- rasterize(v, r)

#plot(x)
#lines(v)
}

\keyword{spatial}

