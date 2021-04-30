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
\S4method{convhull}{SpatVector}(x)
}

\arguments{
  \item{x}{SpatVector}
}

\value{
SpatVector
}

\examples{
p <- vect(system.file("ex/lux.shp", package="terra"))
h <- convhull(p)
}

\keyword{methods}
\keyword{spatial}

