\name{convHull}

\docType{methods}

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
\S4method{convHull}{SpatVector}(x, by="")

\S4method{minRect}{SpatVector}(x, by="")

\S4method{minCircle}{SpatVector}(x, by="")
}

\arguments{
  \item{x}{SpatVector}
  \item{by}{character (variable name), to get a new geometry for groups of input geometries}
}

\value{
SpatVector
}

\examples{
p <- vect(system.file("ex/lux.shp", package="terra"))
h <- convHull(p)

hh <- convHull(p, "NAME_1")
rr <- minRect(p, "NAME_1")

plot(rr, lwd=5, border="gray")
plot(hh, "NAME_1", col=rainbow(10, alpha=.5), lwd=3, add=TRUE, plg=list(x="topright"))
lines(aggregate(p, "NAME_1"), col="blue", lty=2, lwd=2)

}

\keyword{methods}
\keyword{spatial}

