\name{values}

\docType{methods}

\alias{values}
\alias{values<-}
\alias{values,SpatRaster-method}
\alias{values<-,SpatRaster,ANY-method}

\alias{setValues}
\alias{setValues,SpatRaster-method}
\alias{setValues,SpatRaster,ANY-method}

\alias{values,SpatVector-method}
\alias{values<-,SpatVector,data.frame-method}
\alias{values<-,SpatVector,NULL-method}

  
\title{Get or set cell values}

\description{
Get are set all cell values of a SpatRaster or SpatVector. For large SpatRaster objects use \code{\link(init)} instead to set values. 
}

\usage{
\S4method{values}{SpatRaster}(x, mat=TRUE, dataframe=FALSE, row=1, nrows=nrow(x), col=1, ncols=ncol(x))

\S4method{values}{SpatRaster,ANY}(x)<-value

\S4method{setValues}{SpatRaster,ANY}(x, values)

\S4method{values}{SpatVector}(x)

\S4method{values}{SpatVector,data.frame}(x)<-value
}


\arguments{
 \item{x}{SpatRaster or SpatVector}
 \item{mat}{logical. If \code{TRUE}, values are returned as a matrix instead of as a vector, except when dataframe is \code{TRUE}}
 \item{dataframe}{logical. If \code{TRUE}, values are returned as a \code{data.frame} instead of as a vector (also if matrix is \code{TRUE})}

 \item{row}{positive integer. Row number to start from, should be between 1 and nrow(x)}
 \item{nrows}{positive integer. How many rows?}
 \item{col}{positive integer. Column number to start from, should be between 1 and ncol(x)}
 \item{ncols}{positive integer. How many columns? Default is the number of columns left after the start column}


 \item{value}{For SpatRaster: matrix or numeric, the length must match the total number of cells (ncell(x) * nlyr(x)), or be a single value. For SpatVector: data.frame or NULL}
 \item{values}{Same as for \code{value}}
}


\details{
If \code{x} is a \code{SpatRaster}:

If \code{matrix=TRUE}, a matrix is returned in which the values of each layer are represented by a column (with \code{ncell(x)} rows). The values per layer are in cell-order, that is, from top-left, to top-right and then down by row. Use \code{\link{as.matrix}} for an alternative matrix representation where the number of rows and columns matches that of \code{x}, if \code{x} has a single layer.
If \code{matrix=FALSE}, the values are returned as a vector. In cell-order by layer. 

If \code{x} is a \code{SpatVector}: a \code{data.frame}

}


\value{
matrix, vector, SpatRaster, or nothing
}


\seeAlso(\code{\link{init}, \code{\link{focalValues}})

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
