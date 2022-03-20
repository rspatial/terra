\name{unique}

\docType{methods}

\alias{unique}
\alias{unique,SpatRaster-method}
\alias{unique,SpatRaster,ANY-method}
\alias{unique,SpatVector-method}
\alias{unique,SpatVector,ANY-method}

\title{Unique values}

\description{
This function returns the unique values in a SpatRaster or removes duplicates in a SpatVector. 
}

\usage{
\S4method{unique}{SpatRaster}(x, incomparables=FALSE, na.rm=FALSE) 

\S4method{unique}{SpatVector}(x, incomparables=FALSE, ...) 
}

\arguments{
  \item{x}{SpatRaster or SpatVector}
  \item{incomparables}{logical. If \code{FALSE} and \code{x} is a SpatRaster: the unique values are determined for all layers together, and the result is a matrix. If \code{TRUE}, each layer is evaluated separately, and a list is returned. If \code{x} is a SpatVector this argument is as for a \code{data.frame}}
  \item{na.rm}{logical. If \code{TRUE}, \code{NaN} is included if there are any missing values}  
  \item{...}{additional arguments passed on to \code{\link[base]{unique}}}  
}


\value{
If \code{x} is a SpatRaster: data.frame or list (if \code{incomparables=FALSE})

If \code{x} is a SpatVector: SpatVector
}


\examples{
r <- rast(ncols=5, nrows=5)
values(r) <- rep(1:5, each=5)
unique(r)
s <- c(r, round(r/3))
unique(s)
unique(s,TRUE)

v <- vect(cbind(x=c(1:5,1:5), y=c(5:1,5:1)), 
		crs="+proj=utm +zone=1 +datum=WGS84")
nrow(v)
u <- unique(v)
nrow(u)

values(v) <- c(1:5, 1:3, 5:4)
unique(v)
}

\keyword{spatial}

