\name{simplify}

\docType{methods}

\alias{simplify}

\alias{simplify,SpatVector-method}

\title{ 
Simplify geometries
}

\description{
Reduce the number of nodes used to represent geometries.
}

\usage{
\S4method{simplify}{SpatVector}(x, tolerance=0.1)
}

\arguments{
  \item{x}{SpatVector of lines or polygons}
  \item{tolerance}{numeric. The minimum distance between nodes in units of the crs (i.e. degrees for long/lat}
} 

\value{
SpatVector
}


\seealso{\code{\link{sharedPaths}}, \code{\link{gaps}}}


\examples{
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
#w <- simplify(v, .02)
#e <- erase(w)
#g <- gaps(e)
#plot(e, lwd=5, border="light gray")
#polys(g, col="red", border="red")
}


\keyword{methods}
\keyword{spatial}

