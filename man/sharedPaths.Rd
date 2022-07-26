\name{sharedPaths}

\docType{methods}

\alias{sharedPaths}

\alias{sharedPaths,SpatVector-method}

\title{ 
Shared paths
}

\description{
Get shared paths of line or polygon geometries. This can for geometries in a single SpatVector, or between two SpatVectors
}

\usage{
\S4method{sharedPaths}{SpatVector}(x, y=NULL)
}

\arguments{
  \item{x}{SpatVector of lines or polygons}
  \item{y}{missing or SpatVector of lines or polygons}
} 

\value{
SpatVector
}

\seealso{
\code{\link{gaps}}, \code{\link{topology}}
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

z <- sharedPaths(v[3,], v[12,])

}


\keyword{methods}
\keyword{spatial}

