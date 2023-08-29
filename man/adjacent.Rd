\name{adjacent}

\docType{methods}

\alias{adjacent}
\alias{adjacent,SpatRaster-method}
\alias{adjacent,SpatVector-method}

\title{Adjacent cells}

\description{
Identify cells that are adjacent to a set of raster cells. Or identify adjacent polygons
}

\usage{
\S4method{adjacent}{SpatRaster}(x, cells, directions="rook", pairs=FALSE, include=FALSE, symmetrical=FALSE)

\S4method{adjacent}{SpatVector}(x, type="rook", pairs=TRUE, symmetrical=FALSE)
}

\arguments{
  \item{x}{SpatRaster}
  \item{cells}{vector of cell numbers for which adjacent cells should be found. Cell numbers start with 1 in the upper-left corner and increase from left to right and from top to bottom}
  \item{directions}{character or matrix to indicated the directions in which cells are considered connected. The following character values are allowed: "rook" or "4" for the horizontal and vertical neighbors; "bishop" to get the diagonal neighbors; "queen" or "8" to get the vertical, horizontal and diagonal neighbors; or "16" for knight and one-cell queen move neighbors. If \code{directions} is a matrix it should have odd dimensions and have logical (or 0, 1) values}
  \item{pairs}{logical. If \code{TRUE}, a two-column matrix of pairs of adjacent cells is returned. If \code{x} is a \code{SpatRaster} and \code{pairs} is \code{FALSE}, an \code{n*m} matrix is returned where the number of rows \code{n} is \code{length(cells)} and the number of columns \code{m} is the number of neighbors requested with \code{directions}}
  \item{include}{logical. Should the focal cells be included in the result?}
  \item{type}{character. One of "rook", "queen", "touches", or "intersects". "queen" and "touches" are synonyms. "rook" exclude polygons that touch at a single node only. "intersects" includes polygons that touch or overlap}
  \item{symmetrical}{logical. If \code{TRUE} and \code{pairs=TRUE}, an adjacent pair is only included once. For example, if polygon 1 is adjacent to polygon 3, the implied adjacency between 3 and 1 is not reported}
}

\note{
When using global lon/lat rasters, adjacent cells at the other side of the date-line are included.
}


\seealso{\code{\link{relate}}, \code{\link{nearby}}}

\value{
matrix
}

\examples{
r <- rast(nrows=10, ncols=10)
adjacent(r, cells=c(1, 5, 55), directions="queen") 
r <- rast(nrows=10, ncols=10, crs="+proj=utm +zone=1 +datum=WGS84")
adjacent(r, cells=11, directions="rook") 

#same as 
rk <- matrix(c(0,1,0,1,0,1,0,1,0), 3, 3)
adjacent(r, cells=11, directions=rk) 

## note that with global lat/lon data the E and W connect
r <- rast(nrows=10, ncols=10, crs="+proj=longlat +datum=WGS84")
adjacent(r, cells=11, directions="rook") 

f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
a <- adjacent(v, symmetrical=TRUE)
head(a)
}

\keyword{spatial}
