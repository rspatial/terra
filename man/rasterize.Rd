\name{rasterize}

\docType{methods}

\alias{rasterize}
\alias{rasterize,SpatVector,SpatRaster-method}
\alias{rasterize,sf,SpatRaster-method}
\alias{rasterize,matrix,SpatRaster-method}

\title{Rasterize vector data}

\description{
Transfer values associated with the geometries of vector data to a raster
}

\usage{
\S4method{rasterize}{SpatVector,SpatRaster}(x, y, field="", fun, ..., background=NA, touches=FALSE, update=FALSE, 
	cover=FALSE, by=NULL, filename="", overwrite=FALSE, wopt=list())

\S4method{rasterize}{matrix,SpatRaster}(x, y, values=1, fun, ..., background=NA, update=FALSE, 
	by=NULL, filename="", overwrite=FALSE, wopt=list())
}

\arguments{
  \item{x}{SpatVector or a two-column matrix (point coordinates)}
  
  \item{y}{SpatRaster}
  
  \item{field}{character or numeric. If \code{field} is a character, it should a variable name in \code{x}. If \code{field} is numeric it typically is a single number or a vector of length \code{nrow(x)}. The values are recycled to \code{nrow(x)}}

  \item{values}{typically a numeric vector of length \code{1} or \code{nrow(x)}. If the length is below \code{nrow(x)} the values will be recycled to \code{nrow(x)}. Only used when \code{x} is a matrix. Can also be a matrix or data.frame}
  
  \item{fun}{summarizing function for when there are multiple geometries in one cell. For lines and polygons you can only use  \code{"min"}, \code{"max"}, \code{"mean"}, \code{"count"} and \code{"sum"}
  For points you can use any function that returns a single number; for example \code{mean}, \code{length} (to get a count), \code{min} or \code{max}}
  
  \item{...}{additional arguments passed to \code{fun}}
  
  \item{background}{numeric. Value to put in the cells that are not covered by any of the features of \code{x}. Default is \code{NA}}
  
  \item{touches}{logical. If \code{TRUE}, all cells touched by lines or polygons are affected, not just those on the line render path, or whose center point is within the polygon. If \code{touches=TRUE}, \code{add} cannot be \code{TRUE}}

  \item{update}{logical. If \code{TRUE}, the values of the input SpatRaster are updated}
  
  \item{cover}{logical. If \code{TRUE} and the geometry of \code{x} is polygons, the fraction of a cell that is covered by the polygons is returned. This is estimated by determining presence/absence of the polygon in at least 100 sub-cells (more of there are very few cells)} 

  \item{by}{character or numeric value(s) to split \code{x} into multiple groups. There will be a separate layer for each group returned. If \code{x} is a SpatVector, \code{by} can be a column number or name. If \code{x} is a matrix, \code{by} should be a vector that identifies group membership for each row in \code{x}}

  \item{filename}{character. Output filename}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}  
  \item{wopt}{list with additional arguments for writing files as in \code{\link{writeRaster}}}
}

\seealso{
\code{\link{rasterizeGeom}}, \code{\link{rasterizeWin}}, \code{\link{mask}}
}


\value{
SpatRaster
}

\examples{
r <- rast(xmin=0, ncols=18, nrows=18)

# generate points
set.seed(1)
p <- spatSample(r, 1000, xy=TRUE, replace=TRUE)

# rasterize points as a matrix
x <- rasterize(p, r, fun=sum)
y <- rasterize(p, r, value=1:nrow(p), fun=max)

# rasterize points as a SpatVector
pv <- vect(p)
xv <- rasterize(pv, r, fun=sum)

# Polygons
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
r <- rast(v, ncols=75, nrows=100)
z <- rasterize(v, r, "NAME_2")
plot(z)
lines(v)
}

\keyword{spatial}

