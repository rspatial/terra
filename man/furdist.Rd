\name{furdist}

\alias{furdist}
\alias{furdist,SpatVector,SpatVector-method}

\title{Distance to the furthest location on a geometry}

\description{
For each point in \code{y}, compute the distance to the \emph{furthest} location on the edge of the geometries in \code{x}. This is the counterpart of \code{\link{distance}}, which returns the distance to the \emph{nearest} location.

The "edge" is the boundary of polygons, the line itself for lines, and the points for a points geometry.
}

\usage{
\S4method{furdist}{SpatVector,SpatVector}(x, y, pairwise=FALSE, unit="m")
}

\arguments{
  \item{x}{SpatVector of points, lines, or polygons. The geometries whose furthest edge location is measured}
  \item{y}{SpatVector of points. The location(s) from which the distance is measured}
  \item{pairwise}{logical. If \code{TRUE} and if \code{x} and \code{y} have the same number of geometries (or one of them has a single geometry), the pairwise distances (geometry \code{i} of \code{x} with point \code{i} of \code{y}) are returned. If \code{FALSE}, the distances for all combinations are returned}
  \item{unit}{character. Can be either "m" or "km"}
}

\value{
If \code{pairwise=TRUE} a numeric vector for each pair. Otherwise a from-to-distance data.frame. 
}

\seealso{\code{\link{distance}}, \code{\link{nearest}}, \code{\link{snapTo}}}

\examples{
# planar
crs <- "+proj=utm +zone=1"
p <- vect("POLYGON ((0 0, 8 0, 8 9, 0 9, 0 0))", crs=crs)
pts <- vect(cbind(c(1, 5), c(1, 5)), crs=crs)

furdist(p, pts)
furdist(p, pts, unit="km")

# longitude/latitude
pp <- vect("POLYGON ((0 0, 10 0, 10 10, 0 10, 0 0))", crs="lonlat")
q <- vect(cbind(5, 5), crs="lonlat")
furdist(pp, q, unit="km")
}

\keyword{spatial}
