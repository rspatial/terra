\name{split}

\docType{methods}

\alias{split}
\alias{split,SpatRaster,ANY-method}
\alias{split,SpatVector,ANY-method}
\alias{split,SpatVector,SpatVector-method}


\title{Split}

\description{ 
Split a SpatRaster by layer, or a SpatVector by attributes. You can also split the geometry of a polygon SpatVector with another SpatVector. 
}

\usage{
\S4method{split}{SpatRaster,ANY}(x, f)

\S4method{split}{SpatVector,ANY}(x, f)

\S4method{split}{SpatVector,SpatVector}(x, f)
}

\arguments{
  \item{x}{SpatRaster or SpatVector}
  \item{f}{If \code{x} is a SpatRaster: a vector of the length \code{nlyr(x)}. If \code{x} is a SpatVector: a field (variable) name or a vector of the same length as \code{x}; or, if \code{x} is a SpatVector of polygons, a SpatVector of lines or polygons to split the polygon geometries}
}

\value{
list or SpatVector
}


\examples{

## split layers
s <- rast(system.file("ex/logo.tif", package="terra"))   
y <- split(s, c(1,2,1))
sds(y)

## split attributes
v <- vect(system.file("ex/lux.shp", package="terra"))
x <- split(v, "NAME_1")

## split geometries
v <- v[1:5,]
line <- vect(matrix(c(5.79, 6.22, 5.75, 6.1, 5.8, 
	50.14, 50.05, 49.88, 49.85, 49.71), ncol=2), "line")
s <- split(v, line)
}


\keyword{methods}
\keyword{spatial}
