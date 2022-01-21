\name{set_names}

\alias{set_names}
\alias{set_longnames}
\alias{set_varnames}

\alias{set_names,SpatRaster-method}
\alias{set_varnames,SpatRaster-method}
\alias{set_longnames,SpatRaster-method}

\alias{set_names,SpatRasterDataset-method}
\alias{set_varnames,SpatRasterDataset-method}
\alias{set_longnames,SpatRasterDataset-method}

\alias{set_names,SpatVector-method}
\alias{set_names,SpatVectorCollection-method}
  
\title{Change names in-place}

\description{
\code{set_names} is the in-place replacement version of \code{\link{names<-}}. That is, names are assigned to an object without making a copy. That is efficient, but if there is a copy of the object that you made by standard assignment (e.g. with \code{y <- x}), that copy is also changed.
}

\usage{
\S4method{set_names}{SpatRaster}(x, value, index=1:nlyr(x), validate=FALSE)
\S4method{set_varnames}{SpatRaster}(x, value, index=1:nlyr(x))
\S4method{set_longnames}{SpatRaster}(x, value, index=1:nlyr(x))

\S4method{set_names}{SpatRasterDataset}(x, value, index=1:length(x), validate=FALSE)

\S4method{set_names}{SpatVector}(x, value, index=1:ncol(x), validate=FALSE)
}

\arguments{
  \item{x}{SpatRaster}
  \item{value}{character}
  \item{index}{numeric, to select the layer to assign a name tor}
  \item{validate}{logical. Make names valid and/or unique?}
}


\examples{
s <- rast(ncols=5, nrows=5, nlyrs=3)
x <- s
names(s)
names(s) <- c("a", "b", "c")
names(s)
names(x)

x <- s
set_names(s, c("e", "f", "g"))
names(s)
names(x)
}

\keyword{spatial}

