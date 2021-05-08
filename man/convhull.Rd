\name{convhull}

\docType{methods}

\alias{convhull}
\alias{convhull,SpatVector-method}
\alias{convexhull}
\alias{convexhull,SpatVector-method}

\title{ 
Convex hull
}

\description{
Get the convex hull of a vector dataset
}

\usage{
\S4method{convhull}{SpatVector}(x, by="")
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
h <- convhull(p)

hh <- convhull(p, "NAME_1")
plot(hh, "NAME_1", col=rainbow(10, alpha=.5), lwd=3)
lines(aggregate(p, "NAME_1"), col="blue", lty=2, lwd=2)
}

\keyword{methods}
\keyword{spatial}

