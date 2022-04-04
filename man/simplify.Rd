\name{simplifyGeom}

\docType{methods}

\alias{simplifyGeom}

\alias{simplifyGeom,SpatVector-method}

\title{ 
simplifyGeom geometries
}

\description{
Reduce the number of nodes used to represent geometries.
}

\usage{
\S4method{simplifyGeom}{SpatVector}(x, tolerance=0.1, preserveTopology=TRUE)
}

\arguments{
  \item{x}{SpatVector of lines or polygons}
  \item{tolerance}{numeric. The minimum distance between nodes in units of the crs (i.e. degrees for long/lat}
  \item{preserveTopology}{logical. If \code{TRUE} the topology of output geometries is preserved (polygons remain valid)}
} 

\value{
SpatVector
}


\seealso{\code{\link{sharedPaths}}, \code{\link{gaps}}}


\examples{
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
w <- simplifyGeom(v, .02)
e <- erase(w)
g <- gaps(e)
plot(e, lwd=5, border="light gray")
polys(g, col="red", border="red")
}


\keyword{methods}
\keyword{spatial}

