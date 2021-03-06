\name{rasterize}

\docType{methods}

\alias{rasterize}
\alias{rasterize,SpatVector,SpatRaster-method}

\title{Rasterize vector data}

\description{
Transfer vector data to a raster
}

\usage{
\S4method{rasterize}{SpatVector,SpatRaster}(x, y, field="", fun, ..., background=NA, touches=FALSE,
	update=FALSE, sum=FALSE, cover=FALSE, filename="", overwrite=FALSE, wopt=list())
}

\arguments{
  \item{x}{SpatVector}
  \item{y}{SpatRaster}
  
  \item{field}{character or numeric. If \code{field} is a character, it should a variable name in \code{x}. If \code{field} is numeric it typically is a single number or a vector of length \code{nrow(x)}. The values are recycled to \code{nrow(x)}}
  
  \item{fun}{function, summarizing function that returns a single number; for when there are multiple points in one cell. For example \code{mean}, \code{length} (to get a count), \code{min} or \code{max}. Only used if \code{x} consists of points}
  \item{...}{additional arguments passed to \code{fun} if \code{x} has point geometries}
  
  \item{background}{numeric. Value to put in the cells that are not covered by any of the features of \code{x}. Default is \code{NA}}
  
  \item{touches}{logical. If \code{TRUE}, all cells touched by lines or polygons are affected, not just those on the line render path, or whose center point is within the polygon. If \code{touches=TRUE}, \code{add} cannot be \code{TRUE}}

  \item{update}{logical. If \code{TRUE}, the values of the input SpatRaster are updated}
  
  \item{sum}{logical. If \code{TRUE}, the values of overlapping geometries are summed instead of replaced; and \code{background} is set to zero. Only used if \code{x} does not consists of points} 

  \item{cover}{logical. If \code{TRUE} and the geometry of \code{x} is polygons, the fraction of a cell that is covered by the polygons is returned. This is estimated by determining presence/absence of the polygon in at least 100 sub-cells (more of there are very few cells)} 

  \item{filename}{character. Output filename}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}  
  \item{wopt}{list with additional arguments for writing files as in \code{\link{writeRaster}}}
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
r <- rast(v, ncols=75, nrows=100)
#x <- rasterize(v, r, "NAME_2")

#plot(x)
#lines(v)
}

\keyword{spatial}

