\name{convexhull}

\docType{methods}

\alias{convexhull}
\alias{convexhull,SpatVector-method}

\title{ 
Convex hull
}

\description{
Get the convex hull of a vector dataset
}

\usage{
\S4method{convexhull}{SpatVector}(x)
}

\arguments{
  \item{x}{SpatVector}
}

\value{
SpatVector
}

\examples{
p <- vect(system.file("ex/lux.shp", package="terra"))
h <- convexhull(p)
}

\keyword{methods}
\keyword{spatial}

