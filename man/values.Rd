\name{values}

\docType{methods}

\alias{values}
\alias{values,SpatRaster-method}
\alias{values,SpatVector-method}
  
\title{Cell values and geometry attributes}

\description{
Get the cell values of a SpatRaster or the attributes of a SpatVector
}

\usage{
\S4method{values}{SpatRaster}(x, mat=TRUE, dataframe=FALSE, row=1, nrows=nrow(x), col=1, ncols=ncol(x))

\S4method{values}{SpatVector}(x)
}


\arguments{
 \item{x}{SpatRaster or SpatVector}
 \item{mat}{logical. If \code{TRUE}, values are returned as a matrix instead of as a vector, except when dataframe is \code{TRUE}}
 \item{dataframe}{logical. If \code{TRUE}, values are returned as a \code{data.frame} instead of as a vector (also if matrix is \code{TRUE})}

 \item{row}{positive integer. Row number to start from, should be between 1 and nrow(x)}
 \item{nrows}{positive integer. How many rows?}
 \item{col}{positive integer. Column number to start from, should be between 1 and ncol(x)}
 \item{ncols}{positive integer. How many columns? Default is the number of columns left after the start column}
}


\details{
If \code{x} is a \code{SpatRaster}, and \code{matrix=TRUE}, a matrix is returned in which the values of each layer are represented by a column (with \code{ncell(x)} rows). The values per layer are in cell-order, that is, from top-left, to top-right and then down by row. Use \code{\link{as.matrix}} for an alternative matrix representation where the number of rows and columns matches that of \code{x}, if \code{x} has a single layer.
If \code{matrix=FALSE}, the values are returned as a vector. In cell-order by layer. 
}


\value{
matrix or data.frame
}


\seealso{\code{\link{values<-}}, \code{\link{focalValues}}}

\examples{
f <- system.file("ex/elev.tif", package="terra")
r <- rast(f)
r
x <- values(r)
x[3650:3655, ]
r[3650:3655]


ff <- system.file("ex/lux.shp", package="terra")
v <- vect(ff)
y <- values(v)
head(y)
}

\keyword{spatial}
\keyword{methods}
