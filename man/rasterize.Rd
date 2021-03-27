\name{rasterize}

\docType{methods}

\alias{rasterize}
\alias{rasterize,SpatVector,SpatRaster-method}

\title{Rasterize vector data}

\description{
Transfer vector data to a raster
}

\usage{
\S4method{rasterize}{SpatVector,SpatRaster}(x, y, field, values=1, fun, background=NA, touches=FALSE,
	sum=FALSE, cover=FALSE, filename="", ...)
}

\arguments{
  \item{x}{SpatVector}
  \item{y}{SpatRaster}
  
  \item{field}{character or numeric. If \code{field} is a character, it should a variable name in \code{x}. If \code{field} is numeric it should be a single number (index the variable)}

  \item{values}{numeric. Only used if \code{field} is \code{""}. Typically a single number or a vector of length \code{nrow(x)}. The values are recycled to \code{nrow(x)}}
  
  \item{fun}{function, summarizing function that returns a single number; for when there are multiple points in one cell. For example \code{mean}, \code{length} (to get a count), \code{min} or \code{max}. Only used if \code{x} consists of points}
  
  \item{background}{numeric. Value to put in the cells that are not covered by any of the features of \code{x}. Default is \code{NA}}
  
  \item{touches}{logical. If \code{TRUE}, all cells touched by lines or polygons are affected, not just those on the line render path, or whose center point is within the polygon. If \code{touches=TRUE}, \code{add} cannot be \code{TRUE}}
  
  \item{sum}{logical. If \code{TRUE}, the values of overlapping geometries are summed instead of replaced; and \code{background} is set to zero. Only used if \code{x} does not consists of points} 

  \item{cover}{logical. If \code{TRUE} and the geometry of \code{x} is polygons, the fraction of a cell that is covered by the polygons is returned. This is estimated by determining presence/absence of the polygon in at least 100 sub-cells (more of there are very few cells)} 

  \item{filename}{character. Output filename}
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
}

\seealso{
\code{\link{mask}}
}

\note{
To update existing raster data use the output of \code{rasterize} as second argument in \code{\link{cover}} 
}

\value{
SpatRaster
}

\examples{
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
r <- rast(v, ncol=75, nrow=100)
#x <- rasterize(v, r, "NAME_2")

#plot(x)
#lines(v)
}

\keyword{spatial}

