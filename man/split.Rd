\name{split}

\docType{methods}

\alias{split}
\alias{split,SpatRaster-method}
\alias{split,SpatVector-method}


\title{Split}

\description{ 
Split a SpatVector or SpatRaster}

\usage{
\S4method{split}{SpatRaster}(x, f)

\S4method{split}{SpatVector}(x, f)
}

\arguments{
  \item{x}{SpatRaster or SpatVector}
  \item{f}{If \code{x} is a SpatVector: a field (variable) name or a vector of the same length as \code{x}. If \code{x} is a SpatRaster: a vector of the length \code{nlyr(x)}}
}

\value{
Same as \code{x}
}


\examples{
v <- vect(system.file("ex/lux.shp", package="terra"))
x <- split(v, "NAME_1")

s <- rast(system.file("ex/logo.tif", package="terra"))   
y <- split(s, c(1,2,1))
sds(y)
}


\keyword{methods}
\keyword{spatial}
