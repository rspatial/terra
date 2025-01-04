\name{convHull}

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
Convex hull, minimal bounding rotated rectangle, and minimal bounding circle
}

\description{
Get the convex hull, the minimal bounding rotated rectangle, or minimal bounding circle of a SpatVector
}

\usage{
\S4method{hull}{SpatVector}(x, type="convex", by="")
}

\arguments{
  \item{x}{SpatVector}
  \item{type}{character. One of "convex", "rectangle", "minCircle"}
  \item{by}{character (variable name), to get a new geometry for groups of input geometries}
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

