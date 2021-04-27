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
\S4method{adjacent}{SpatRaster}(x, cells, directions="rook", include=FALSE)

\S4method{adjacent}{SpatVector}(x, type="rook", pairs=TRUE, symmetrical=FALSE	)
}

\arguments{
  \item{x}{SpatRaster}
  \item{cells}{vector of cell numbers for which adjacent cells should be found. Cell numbers start with 1 in the upper-left corner and increase from left to right and from top to bottom}
  \item{directions}{the directions in which cells should be connected: "rook" (4 directions), "queen" (8 directions), "16" (knight and one-cell queen moves), or "bishop" to connect cells with one-cell diagonal moves. }
  \item{include}{logical. Should the focal cells be included in the result? }
  \item{type}{character. One of "rook", "queen", "touches", or "intersects". "queen" and "touches" are synonyms. "rook" exclude polygons that touch at a single node only. "intersects" includes polygons that touch or overlap}
  \item{pairs}{logical. If \code{TRUE} a "from", "to" matrix is returned}   
  \item{symmetrical}{logical. If \code{TRUE}, an adjacent pair is only included once. For example, if polygon 1 is adjacent to polygon 3, the implied adjacency between 3 and 1 is not reported}
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
# global lat/lon wraps around
r <- rast(nrows=10, ncols=10, crs="+proj=longlat +datum=WGS84")
adjacent(r, cells=11, directions="rook") 

f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
a <- adjacent(v, symmetrical=TRUE)
head(a)
}

\keyword{spatial}
