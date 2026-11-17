\name{hull}

\docType{methods}

\alias{hull}
\alias{hull,SpatVector-method}
\alias{convHull}
\alias{convHull,SpatVector-method}

\title{ 
Convex, concave, rectangular and circular hulls
}

\description{
Compute hulls around SpatVector geometries. This can be the convex hull, the minimal bounding rotated rectangle, the minimal bounding circle, or a concave hull. The concaveness of the concave hull can be specified in different ways.

The old method \code{convHull} is deprecated and will be removed in a future version. 
}

\usage{
\S4method{hull}{SpatVector}(x, type="convex", by="", param=1, allowHoles=TRUE, tight=TRUE, percent=100)
}

\arguments{
  \item{x}{SpatVector}
  \item{type}{character. One of "convex", "rectangle", "circle", "concave_ratio", "concave_length"}
  \item{by}{character (variable name), to get a new geometry for groups of input geometries. Use \code{""} to get a hull around all geometries combined. Use \code{NA} to get a hull around each geometry in \code{x}}
  \item{param}{numeric. For \code{type="concave_ratio"}, use an edge length ratio between 0 and 1. For \code{type="concave_length"}, the meaning depends on the input geometry: use a ratio between 0 and 1 for polygons, or a maximum edge length in coordinate units (0 or greater) for points and lines}
  \item{allowHoles}{logical. May the output polygons contain holes? For "concave_*" methods only}
  \item{tight}{logical. Should the hull follow the outer boundaries of the input polygons? For "concave_length" with polygon geometry only}
  \item{percent}{numeric > 0 and <= 100. If below 100, only that percentage of geometries closest to the centroid is used}
}

\value{
SpatVector
}


\details{
A concave hull is a polygon which contains all the points of the input. It can be a better representation of the input data (typically points) than the convex hull. There are many possible concave hulls with different degrees of concaveness. These can be created with argument \code{param}.

For \code{type="concave_ratio"}, larger ratios retain longer edges. A ratio of 0 gives the most concave hull and 1 gives the convex hull. This method requires GEOS 3.11 or later.

For \code{type="concave_length"} with polygons, the ratio applies to the edges between polygons. At 0 the original polygons are retained. The \code{tight} argument controls whether the hull follows their outer boundaries. This method requires GEOS 3.11 or later. For points and lines, a larger maximum edge length gives a less concave hull. This method requires GEOS 3.12 or later.
}



\examples{
p <- vect(system.file("ex/lux.shp", package="terra"))
h <- hull(p)

plot(p)
lines(h, col="orange")

hh <- hull(p, "convex", by="NAME_1")
lines(hh, col="purple")

pts <- centroids(p)
plot(pts, ext=ext(p)+0.1)
lines(hull(pts, type="convex"), col="darkgreen")
lines(hull(pts, type="rectangle"), col="blue")
lines(hull(pts, type="circle"), col="red")
lines(hull(pts, percent=60), col="orange")

}

\keyword{methods}
\keyword{spatial}

