\name{unique}

\docType{methods}

\alias{unique}
\alias{unique,SpatRaster,ANY-method}

\title{Unique values}

\description{
This function returns the unique values in a SpatRaster.
}

\usage{
\S4method{unique}{SpatRaster}(x, incomparables=FALSE, ...) 
}

\arguments{
  \item{x}{SpatRaster}
  \item{incomparables}{logical. If \code{TRUE}, the unique values are determined for all layers together, and the result is a matrix. If \code{FALSE}, each layer is evaluated separately, and a list is returned.}
  \item{...}{additional arguments. none implemented}  
}


\value{
vector or matrix
}


\examples{
r <- rast(ncol=5, nrow=5)
values(r) <- rep(1:5, each=5)
unique(r)
s <- c(r, round(r/3))
unique(s)
}

\keyword{spatial}

