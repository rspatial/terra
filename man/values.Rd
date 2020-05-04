\name{values}

\docType{methods}


\alias{values}
\alias{values<-}
\alias{values,SpatRaster-method}
\alias{values<-,SpatRaster,ANY-method}
\alias{values,SpatVector-method}
\alias{values<-,SpatVector,data.frame-method}
  
\title{Get cell values}

\description{
values returns all cell values of a SpatRaster (a matrix), or all the attributes of a SpatVector (a data.frame).
}

\usage{
\S4method{values}{SpatRaster}(x, mat=TRUE, ...)
\S4method{values}{SpatRaster,ANY}(x, ...)<-value
\S4method{values}{SpatVector}(x, ...)
\S4method{values}{SpatVector,data.frame}(x, ...)<-value
}

\arguments{
 \item{x}{SpatRaster or SpatVector}
 \item{mat}{logical. If \code{TRUE}, values are returned as a matrix instead of as a vector}
 \item{value}{For SpatRaster: matrix or numeric, the length must match the total number of cells (ncell(x) * nlyr(x)), or be a single value. For SpatVector: data.frame}
 \item{...}{additional arguments. none implemented}
}


\details{
If \code{x} is a \code{SpatRaster}:

If \code{matrix=TRUE}, a matrix is returned in which the values of each layer are represented by a column (with \code{ncell(x)} rows). The values per layer are in cell-order, that is, from top-left, to top-right and then down by row. Use \code{\link{as.matrix}} for an alternative matrix representation where the number of rows and columns matches that of \code{x}, if \code{x} has a single layer.
If \code{matrix=FALSE}, the values are returned as a vector. In cell-order by layer. 

If \code{x} is a \code{SpatVector}: a \code{data.frame}

}


\value{
matrix or vector
}


\examples{
r <- rast(system.file("ex/test.tif", package="terra"))
r
v <- values(r)
values(r) <- v * 10
r

f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
x <- values(v)
x
values(v) <- x[,1:2]
v
}

\keyword{spatial}
\keyword{methods}
