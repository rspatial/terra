\name{snapTo}

\alias{snapTo}
\alias{snapTo,SpatVector,SpatVector-method}

\title{Snap points to the edge of lines or polygons}

\description{
Move the points in \code{x} onto the edge of the corresponding line or polygon in \code{y}.
A point that touches the the target geometry is not changed. A point that is outside the geometry is moved to the nearest location on the edge (the boundary of a polygon, or the line itself). 
}

\usage{
\S4method{snapTo}{SpatVector,SpatVector}(x, y, paired=TRUE, ...)
}

\arguments{
  \item{x}{SpatVector of points}
  \item{y}{SpatVector of lines or polygons. Either with the same number of geometries as there are points in \code{x}, or with a single geometry}
  \item{paired}{logical. If \code{TRUE} the points in \code{x} are paired with the geometries in \code{y} (point \code{i} with geometry \code{i}). \code{y} may instead have a single geometry, in which case all points are snapped to that geometry. Otherwise the points in \code{x} are snapped to any geometry in \code{y}}
  \item{...}{additional arguments. None implemented}
}

\value{
SpatVector of points
}

\seealso{\code{\link{nearest}}, \code{\link{distance}}, \code{\link{furdist}}}

\examples{
crs <- "+proj=utm +zone=1"
p <- vect("POLYGON ((0 0, 8 0, 8 9, 0 9, 0 0))", crs=crs)

# one point inside the polygon (unchanged), one outside (snapped to the edge)
pts <- vect(cbind(c(4, 12), c(4, 4)), crs=crs)
s <- snapTo(pts, p)
crds(s)

# snap to a line
L <- vect("LINESTRING (0 0, 10 0)", crs=crs)
q <- vect(cbind(c(5, 5), c(0, 3)), crs=crs)
crds(snapTo(q, L))
}

\keyword{spatial}
