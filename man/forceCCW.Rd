\name{forceCCW}

\alias{forceCCW}
\alias{forceCCW,SpatVector-method}


\title{force counter-clockwise polygons}


\description{
Assure that the nodes of outer rings of polygons are in counter-clockwise order.
}

\usage{
\S4method{forceCCW}{SpatVector}(x)
}

\arguments{
  \item{x}{SpatVector of polygons}
}

\value{
SpatVector
}

\examples{
p <- vect("POLYGON ((2 45, 2 55, 18 55, 18 45, 2 45))")
pcc <- forceCCW(p)
geom(pcc, wkt=TRUE)
}

\keyword{spatial}

