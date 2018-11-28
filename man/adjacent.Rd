\name{adjacent}


\docType{methods}

\alias{adjacent}
\alias{adjacent,SpatRaster-method}

\title{Adjacent cells}

\description{
Identify cells that are adjacent to a set of cells on a raster. 
}

\usage{
\S4method{adjacent}{SpatRaster}(x, cells, directions, include, ...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{cells}{vector of cell numbers for which adjacent cells should be found. Cell numbers start with 1 in the upper-left corner and increase from left to right and from top to bottom}
  \item{directions}{The directions in which cells should be connected: "rook" (4 directions), "queen" (8 directions), 16 (knight and one-cell queen moves), or "bishop" to connect cells with one-cell diagonal moves. }
  \item{include}{logical. Should the focal cells be included in the result? }
}


\value{
vector with adjacent cells. 
}

\examples{
r <- rast(nrows=10, ncols=10)
adjacent(r, cells=c(5, 55), directions="queen") 
}

\keyword{spatial}
