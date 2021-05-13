\name{topology}

\docType{methods}

\alias{sharedPaths}

\alias{sharedPaths,SpatVector-method}

\title{ 
sharedPaths
}

\description{
Get the shared paths of lines or polygons
}

\usage{
\S4method{sharedPaths}{SpatVector}(x)
}

\arguments{
  \item{x}{SpatVector of lines or polygons}
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

