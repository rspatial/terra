\name{unique}

\docType{methods}

\alias{unique}
\alias{unique,SpatRaster-method}
\alias{unique,SpatRaster,ANY-method}
\alias{unique,SpatVector-method}

\title{Unique values}

\description{
This function returns the unique values in a SpatRaster.
}

\usage{
\S4method{unique}{SpatRaster}(x, incomparables=FALSE, ...) 

\S4method{unique}{SpatVector}(x, incomparables=FALSE, ...) 
}

\arguments{
  \item{x}{SpatRaster or SpatVector}
  \item{incomparables}{logical. If \code{TRUE} and \code{x} is a SpatRaster: the unique values are determined for all layers together, and the result is a matrix. If \code{FALSE}, each layer is evaluated separately, and a list is returned. If \code{x} is a SpatVector this argument is as for a data.frame}
  \item{...}{additional arguments passed on to base::unique (SpatVector only)}  
}


\value{
If \code{x} is a SpatRaster: vector or matrix

If \code{x} is a SpatVector: SpatVector
}


\examples{
r <- rast(ncol=5, nrow=5)
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

}

\keyword{spatial}

