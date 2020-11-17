\name{separate}

\docType{methods}

\alias{separate}
\alias{separate,SpatRaster-method}

\title{separate}

\description{
Create a SpatRaster with a layer for each class (value, or subset of the values) in the input SpatRaster. For example, if the input has vegetation types, this function will create a layer (presence/absence; dummy variable) for each of these classes. Classes and cell values are always truncated to integers.
}

\usage{
\S4method{separate}{SpatRaster}(x, classes=NULL, keep=FALSE, other=0, filename="", overwrite=FALSE, wopt=list(),...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{classes}{numeric. The values (classes) for which layers should be made. If \code{NULL} all classes are used}
  \item{keep}{logical. If \code{TRUE}, cells that are of the class represented by a layer get that value, rather than a value of 1}
  \item{other}{numeric. Value to assign to cells that are not of the class represented by a layer}
  \item{filename}{character. Output filename. Optional}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list. Options for writing files as in \code{\link{writeRaster}}}
  \item{...}{additional arguments. None implemented}    
}

\value{
SpatRaster
}

\examples{
r <- rast(nrow=5, ncol=5)
values(r) <- rep(c(1:4, NA), each=5)
b <- separate(r)
bb <- separate(r, keep=TRUE, other=NA)
}

\keyword{spatial}
