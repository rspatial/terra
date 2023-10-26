\name{rotate}

\docType{methods}

\alias{rotate}
\alias{rotate,SpatRaster-method}
\alias{rotate,SpatVector-method}

\title{Rotate data along longitude}

\description{
Rotate a SpatRaster that has longitude coordinates from 0 to 360, to standard coordinates between -180 and 180 degrees (or vice-versa). Longitude between 0 and 360 is frequently used in global climate models.

Rotate a SpatVector as for a SpatRaster \code{split=TRUE}, or to correct for coordinates that are connected across the date line (and end up at the "other side" of the longitude scale) are reconnected.
}

\usage{
\S4method{rotate}{SpatRaster}(x, left=TRUE, filename="", ...)

\S4method{rotate}{SpatVector}(x, longitude=0, split=FALSE, left=TRUE, normalize=FALSE)
}

\arguments{
 \item{x}{SpatRaster or SpatVector}
 \item{left}{logical. If \code{TRUE}, rotate to the left, else to the right}
 \item{filename}{character. Output filename}
 \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}

 \item{longitude}{numeric. The longitude around which to rotate}
 \item{split}{logical. Should geometries be split at \code{longitude}?}
 \item{normalize}{logical. Should the output be normalized to longitudes between -180 and 180? See \code{\link{normalize.longitude}}}
}

\value{
SpatRaster
}

\seealso{\code{\link{shift}} and \code{\link{spin}}}

\examples{
x <- rast(nrows=9, ncols=18, nl=3, xmin=0, xmax=360)
v <- rep(as.vector(t(matrix(1:ncell(x), nrow=9, ncol=18))), 3)
values(x) <- v
z <- rotate(x)

\dontrun{
#SpatVector
p <- rbind(c(3847903, 1983584 ), c(3847903, 5801864), c(8301883, 5801864), c(8301883, 1983584 ))
p <- vect(p, "polygons", crs="+init=EPSG:3347")
d <- densify(p, 100000)
g <- project(d, "+proj=longlat")

x <- rotate(g, 50)
plot(g)
lines(x, col="red")
}

## rotate countries to 0-360 longitude
#w <- geodata::world(path=".")
#x <- rotate(w, long=0, split=TRUE, left=FALSE)

}

\keyword{spatial}
