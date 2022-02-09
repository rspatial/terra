\name{inplace}

\alias{set.ext}
\alias{set.ext,SpatRaster-method}
\alias{set.ext,SpatVector-method}

\alias{set.values}
\alias{set.values,SpatRaster-method}

\alias{set.cats}
\alias{set.cats,SpatRaster-method}

\alias{set.names}
\alias{set.names,SpatRaster-method}
\alias{set.names,SpatRasterDataset-method}
\alias{set.names,SpatVector-method}
\alias{set.names,SpatVectorCollection-method}

\alias{set.crs}
\alias{set.crs,SpatRaster-method}
\alias{set.crs,SpatVector-method}

  
\title{Change values in-place}

\description{
These "in-place" replacement methods assign new value to an object without making a copy. That is efficient, but if there is a copy of the object that you made by standard assignment (e.g. with \code{y <- x}), that copy is also changed.

\code{set.names} is the in-place replacement version of \code{\link{names<-}}. 
\code{set.ext} is the in-place replacement version of \code{\link{ext<-}}
\code{set.values} is the in-place replacement version of \code{\link{[<-}}
\code{set.cats} is the in-place replacement version of \code{\link{categories}}
\code{set.crs} is the in-place replacement version of \code{\link{crs<-}}
}

\usage{
\S4method{set.names}{SpatRaster}(x, value, index=1:nlyr(x), validate=FALSE)
\S4method{set.names}{SpatRasterDataset}(x, value, index=1:length(x), validate=FALSE)
\S4method{set.names}{SpatVector}(x, value, index=1:ncol(x), validate=FALSE)

\S4method{set.ext}{SpatRaster}(x, value)
\S4method{set.ext}{SpatVector}(x, value)

\S4method{set.crs}{SpatRaster}(x, value)
\S4method{set.crs}{SpatVector}(x, value)

\S4method{set.values}{SpatRaster}(x, cells, values)

\S4method{set.cats}{SpatRaster}(x, layer=1, value, index)
}

\arguments{
  \item{x}{SpatRaster}
  \item{value}{character (\code{set.names}). For \code{set.cats}: a data.frame with columns (value, category) or vector with category names}
  \item{index}{positive integer indicating layer(s) to assign a name to, or the index to select the active category}
  \item{validate}{logical. Make names valid and/or unique?}

  \item{cells}{cell numbers}
  \item{values}{replacement values}

  \item{layer}{positive integer indicating to which layer to you want to assign these categories}
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

