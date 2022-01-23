\name{inplace}

\alias{set.ext}
\alias{set.ext,SpatRaster-method}
\alias{set.ext,SpatVector-method}

\alias{set.values}
\alias{set.values,SpatRaster-method}

\alias{set.names}
\alias{set.names,SpatRaster-method}
\alias{set.names,SpatRasterDataset-method}
\alias{set.names,SpatVector-method}
\alias{set.names,SpatVectorCollection-method}
  
\title{Change values in-place}

\description{
These "in-place" replacement methods assign new value to an object without making a copy. That is efficient, but if there is a copy of the object that you made by standard assignment (e.g. with \code{y <- x}), that copy is also changed.

\code{set.names} is the in-place replacement version of \code{\link{names<-}}. 
\code{set.ext} is the in-place replacement version of \code{\link{ext<-}}
\code{set.values} is the in-place replacement version of \code{\link{[<-}}
}

\usage{
\S4method{set.names}{SpatRaster}(x, value, index=1:nlyr(x), validate=FALSE)
\S4method{set.names}{SpatRasterDataset}(x, value, index=1:length(x), validate=FALSE)
\S4method{set.names}{SpatVector}(x, value, index=1:ncol(x), validate=FALSE)

\S4method{set.ext}{SpatRaster}(x, value)
\S4method{set.ext}{SpatVector}(x, value)

\S4method{set.values}{SpatRaster}(x, cells, values)
}

\arguments{
  \item{x}{SpatRaster}
  \item{value}{character}
  \item{index}{numeric, to select the layer to assign a name to}
  \item{validate}{logical. Make names valid and/or unique?}

  \item{cells}{cell numbers}
  \item{values}{replacement values}
}


\examples{
s <- rast(ncols=5, nrows=5, nlyrs=3)
x <- s
names(s)
names(s) <- c("a", "b", "c")
names(s)
names(x)

x <- s
set.names(s, c("e", "f", "g"))
names(s)
names(x)

set.ext(x, c(0,180,0,90))

set.values(x, 1:10, 5)
}

\keyword{spatial}

