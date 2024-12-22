\name{nseg}

\alias{nseg}
\alias{nseg,SpatVector-method}

  
\title{Number of segments}

\description{
Count the number of segements in a SpatVector of lines or polygons
}

\usage{
\S4method{nseg}{SpatVector}(x)
}


\arguments{
  \item{x}{SpatVector}
}

\value{
numeric
}

\examples{
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
nseg(v)
}

\keyword{methods}
\keyword{spatial}
