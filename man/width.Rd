\name{width}

\docType{methods}

\alias{width}
\alias{clearance}

\alias{width,SpatVector-method}
\alias{clearance,SpatVector-method}

\title{ 
SpatVector geometric properties
}

\description{
Geometric properties
}

\usage{
\S4method{width}{SpatVector}(x)
\S4method{clearance}{SpatVector}(x)
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
width(v)
clearance(v)
}


\keyword{methods}
\keyword{spatial}

