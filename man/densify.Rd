\name{densify}

\docType{methods}

\alias{densify}
\alias{densify,SpatVector-method}

\title{ 
Add additional nodes to lines or polygons
}

\description{
Add additional nodes to lines or polygons. This can be useful to do prior to using \code{project} such that the path does not change too much.
}

\usage{
\S4method{densify}{SpatVector}(x, interval, equalize=TRUE)
}

\arguments{
  \item{x}{SpatVector}
  \item{interval}{numeric larger than 1, specifying the desired minimum interval between nodes}
  \item{equalize}{logical. If \code{TRUE}, new nodes are spread at equal intervals between old nodes}
} 

\value{
SpatVector
}


\examples{
v <- vect(rbind(c(-120,-20), c(-80,5), c(-40,-60), c(-120,-20)), type="polygons", crs="+proj=lonlat")
vd <- densify(v, 200000)

p <- project(v, "+proj=robin")
pd <- project(vd, "+proj=robin")

# good 
plot(pd, col="gray", border="red", lwd=10)
points(pd, col="gray")

# bad
lines(p, col="blue", lwd=3)
points(p, col="blue", cex=2)
plot(p, col="blue", alpha=.1, add=TRUE)
legend("topright", c("good", "bad"), col=c("red", "blue"), lty=1, lwd=3)

## the other way around does not work
## unless the original data was truly planar (e.g. derived from a map)
x <- densify(p, 250000)
y <- project(x, "+proj=longlat")
# bad
plot(y)
# good
lines(vd, col="red")
}

\keyword{methods}
\keyword{spatial}

