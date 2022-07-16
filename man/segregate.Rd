\name{segregate}

\docType{methods}

\alias{segregate}
\alias{segregate,SpatRaster-method}

\title{segregate}

\description{
Create a SpatRaster with a layer for each class (value, or subset of the values) in the input SpatRaster. For example, if the input has vegetation types, this function will create a layer (presence/absence; dummy variable) for each of these classes.

This is called "one-hot encoding" or "dummy encoding" (for a dummy encoding scheme you can remove (any) one of the output layers as it is redundant).
}

\usage{
\S4method{segregate}{SpatRaster}(x, classes=NULL, keep=FALSE, other=0, round=FALSE, digits=0, filename="", ...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{classes}{numeric. The values (classes) for which layers should be made. If \code{NULL} all classes are used}
  \item{keep}{logical. If \code{TRUE}, cells that are of the class represented by a layer get that value, rather than a value of 1}
  \item{other}{numeric. Value to assign to cells that are not of the class represented by a layer}
  \item{round}{logical. Should the values be rounded first?}
  \item{digits}{integer. Number of digits to round the values to}
  \item{filename}{character. Output filename}
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
}

\value{
SpatRaster
}

\examples{
r <- rast(nrows=5, ncols=5)
values(r) <- rep(c(1:4, NA), each=5)
b <- segregate(r)
bb <- segregate(r, keep=TRUE, other=NA)
}

\keyword{spatial}
