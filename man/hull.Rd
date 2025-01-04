\name{hull}

\docType{methods}

\alias{hull}
\alias{hull,SpatVector-method}
\alias{convHull}
\alias{convHull,SpatVector-method}
\alias{minRect}
\alias{minRect,SpatVector-method}
\alias{minCircle}
\alias{minCircle,SpatVector-method}

\title{ 
Convex, concave, rectangular and circular hulls
}

\description{
Compute hulls around SpatVector geometries. This can be the convex hull, the minimal bounding rotated rectangle, the minimal bounding circle, or a concave hull. The concaveness of the concave hull can be specified in different ways.
}

\usage{
\S4method{hull}{SpatVector}(x, type="convex", by="", param=1, allowHoles =TRUE, tight = TRUE)
}

\arguments{
  \item{x}{SpatVector}
  \item{type}{character. One of "convex", "rectangle", "circle", "concave_ratio", "concave_length", "concave_polygons"}
  \item{by}{character (variable name), to get a new geometry for groups of input geometries}
  \item{param}{numeric between 0 and 1. For the "concave_*" types only. For \code{type="concave_ratio"} this is The edge length ratio value, between 0 and 1. For \code{type="concave_length"} this the maximum edge length (a value > 0). For \code{type="concave_polygons"} thism specifies the maximum Edge Length as a fraction of the difference between the longest and shortest edge lengths between the polygons. This normalizes the maximum edge length to be scale-free. A value of 1 produces the convex hull; a value of 0 produces the original polygons}
  \item{allowHoles}{logical. May the polygonal output contain holes? For the "concave_* types only}
  \item{tight}{logical. Should the hull follow the outer boundaries of the input polygons? For the "concave_* types only}
}

\value{
SpatVector
}

\examples{
p <- vect(system.file("ex/lux.shp", package="terra"))
h <- hull(p)

hh <- hull(p, "convex", by="NAME_1")

geos_version <- libVersion("geos", TRUE)

if ((geos_version$major >= 3) && (geos_version$minor >= 6)) {
	rr <- hull(p, "rectangle", "NAME_1")
	plot(rr, lwd=5, border="gray")
	plot(hh, "NAME_1", col=rainbow(10, alpha=.5), lwd=3, add=TRUE, plg=list(x="topright"))
	lines(aggregate(p, "NAME_1"), col="blue", lty=2, lwd=2)
}

}

\keyword{methods}
\keyword{spatial}

