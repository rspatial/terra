\name{convexHull}

\docType{methods}

\alias{convexHull}
\alias{convexHull,SpatVector-method}
\alias{minRect}
\alias{minRect,SpatVector-method}

\title{ 
Convex hull and minimal rotated rectangle
}

\description{
Get the convex hull or the minimal rotated rectangle of a SpatVector
}

\usage{
\S4method{convexHull}{SpatVector}(x, by="")

\S4method{minRect}{SpatVector}(x, by="")
}

\arguments{
  \item{x}{SpatVector}
  \item{by}{character (variable name), to make convex hulls by group}
}

\value{
SpatVector
}

\examples{
p <- vect(system.file("ex/lux.shp", package="terra"))
h <- convexHull(p)

hh <- convexHull(p, "NAME_1")
rr <- minRect(p, "NAME_1")

plot(rr, lwd=5, border="gray")
plot(hh, "NAME_1", col=rainbow(10, alpha=.5), lwd=3, add=TRUE, plg=list(x="topright"))
lines(aggregate(p, "NAME_1"), col="blue", lty=2, lwd=2)

}

\keyword{methods}
\keyword{spatial}

