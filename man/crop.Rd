\name{crop}

\docType{methods}

\alias{crop}
\alias{crop,SpatRaster-method}
\alias{crop,SpatRaster,ANY-method}
\alias{crop,SpatVector,ANY-method}
\alias{crop,SpatVector,SpatVector-method}

\title{Cut out a geographic subset}

\description{
Cut out a part of a SpatRaster with a SpatExtent, or another object from which an extent can be obtained. With a SpatRaster you can only extract rectangular areas, but see \code{\link{mask}} for setting cell values within SpatRaster to \code{NA}.

You can crop a SpatVector with a rectangle, or with another vector (if these are not polygons, the minimum convex hull is used)
}

\usage{
\S4method{crop}{SpatRaster}(x, y, snap="near", filename="", ...)

\S4method{crop}{SpatVector,ANY}(x, y)

\S4method{crop}{SpatVector,SpatVector}(x, y)
}

\arguments{
  \item{x}{SpatRaster or SpatVector}
  \item{y}{SpatExtent or other object that has a SpatExtent (\code{\link{ext}} returns a SpatExtent), or a SpatVector to crop another SpatVector}
  \item{snap}{character. One of "near", "in", or "out"}  
  \item{filename}{character. Output filename}
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
}


\value{
SpatRaster
}


\examples{
r <- rast(xmin=0, xmax=10, ymin=0, ymax=10, nrows=25, ncols=25)
values(r) <- 1:ncell(r)
e <- ext(-5, 5, -5, 5)
rc <- crop(r, e)

# vector
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
e <- ext(6.15, 6.3, 49.7, 49.8)
x <- crop(v, e)
#plot(x, "NAME_1")
}

\keyword{spatial}

