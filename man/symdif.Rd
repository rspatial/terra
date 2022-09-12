\name{symdif}

\docType{methods}

\alias{symdif}
\alias{symdif,SpatVector,SpatVector-method}

\title{ 
Symmetrical difference
}

\description{
Symmetrical difference of polygons
}
\usage{
\S4method{symdif}{SpatVector,SpatVector}(x, y)
}

\arguments{
  \item{x}{SpatVector}
  \item{y}{SpatVector}
}

\value{
SpatVector
}


\seealso{
\code{\link{erase}} 
}

\examples{
p <- vect(system.file("ex/lux.shp", package="terra"))
b <- as.polygons(ext(6, 6.4, 49.75, 50))
#sd <- symdif(p, b)
#plot(sd, col=rainbow(12))
}

\keyword{methods}
\keyword{spatial}

