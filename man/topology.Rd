\name{topology}

\docType{methods}

\alias{sharedPaths}
\alias{mergeLines}
\alias{snap}
\alias{makeNodes}
\alias{removeDupNodes}

\alias{sharedPaths,SpatVector-method}
\alias{mergeLines,SpatVector-method}
\alias{snap,SpatVector-method}
\alias{makeNodes,SpatVector-method}
\alias{removeDupNodes,SpatVector-method}

\title{ 
Vector topology methods
}

\description{
Get the shared paths of lines or polygons
}

\usage{
\S4method{sharedPaths}{SpatVector}(x)
\S4method{snap}{SpatVector}(x, tolerance)
\S4method{mergeLines}{SpatVector}(x)
\S4method{removeDupNodes}{SpatVector}(x, digits = -1)
\S4method{makeNodes}{SpatVector}(x)
}

\arguments{
  \item{x}{SpatVector of lines or polygons}
  \item{tolerance}{numeric}
  \item{digits}{numeric}
} 

\value{
SpatVector
}


\examples{
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
plot(v, col="light gray")
text(v, halo=TRUE)

x <- sharedPaths(v)
lines(x, col="red", lwd=2)
text(x, col="blue", halo=TRUE, cex=0.8)
head(x)
}


\keyword{methods}
\keyword{spatial}

