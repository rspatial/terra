\name{voronoi}

\alias{delaunay}
\alias{voronoi}
\alias{voronoi,SpatVector-method}
\alias{delaunay,SpatVector-method}


\title{Voronoi diagram and Delaunay triangles}

\description{
Get a Voronoi diagram or Delaunay triangles for points, or the nodes of lines or polygons
}

\usage{
\S4method{voronoi}{SpatVector}(x, bnd=NULL, tolerance=0, as.lines=FALSE, deldir=FALSE) 

\S4method{delaunay}{SpatVector}(x, tolerance=0, as.lines=FALSE) 
}


\arguments{
\item{x}{SpatVector}
\item{bnd}{SpatVector to set the outer boundary of the voronoi diagram}
\item{tolerance}{numeric >= 0, snapping tolerance (0 is no snapping)}
\item{as.lines}{logical. If \code{TRUE}, lines are returned without the outer boundary}
\item{deldir}{logical. If \code{TRUE}, the \code{\link[deldir]{deldir}} is used instead of the GEOS C++ library method. It has been reported that \code{deldir} does not choke on very large data sets}
}


\value{SpatVector}


\examples{
wkt <- c("MULTIPOLYGON ( ((40 40, 20 45, 45 30, 40 40)), 
  ((20 35, 10 30, 10 10, 30 5, 45 20, 20 35),(30 20, 20 15, 20 25, 30 20)))", 
  "POLYGON ((0 -5, 10 0, 10 -10, 0 -5))")
x <- vect(wkt)
v <- voronoi(x)
v

d <- delaunay(x)
d

plot(v, lwd=2, col=rainbow(15))
lines(x, col="gray", lwd=2)
points(x)
}

\keyword{spatial}

