\name{values}

\docType{methods}


\alias{values}
\alias{values,SpatRaster-method}
\alias{values<-}
\alias{values<-,SpatRaster-method}
  
\title{Get cell values}

\description{
values returns all cell values of a SpatRaster object.
}

\usage{
\S4method{values}{SpatRaster}(x, matrix=TRUE, ...)
\S4method{values}{SpatRaster}(x,...)<-value
}

\arguments{
 \item{x}{SpatRaster object}
 \item{matrix}{logical. If \code{TRUE}, values are returned as a matrix instead of as a vector}
 \item{value}{matrix or numeric, the length must match the total number of cells (ncell(x) * nlyr(x)), or be a single value}
 \item{...}{additional arguments. none implemented}
}


\details{
If \code{matrix=TRUE}, a matrix is returned in which the values of each layer are represented by a column (with \code{ncell(x)} rows). The values per layer are in cell-order, that is, from top-left, to top-right and then down by row. Use \code{\link{as.matrix}} for an alternative matrix representation where the number of rows and columns matches that of \code{x}, if \code{x} has a single layer.
If \code{matrix=FALSE}, the values are returned as a vector. In cell-order by layer. 
}


\value{
matrix or vector
}


\examples{
r <- rast(system.file("exdata/test.tif", package="terra"))
r
v <- values(r)
values(r) <- v * 10
r
}

\keyword{spatial}
\keyword{methods}
