\name{gaps}

\docType{methods}

\alias{gaps}
\alias{gaps,SpatVector-method}
\alias{gaps,SpatVector,SpatExtent-method}

\title{Find gaps between polygons}

\description{
Get the gaps between polygons of a SpatVector
}

\usage{
\S4method{gaps}{SpatVector}(x)
}

\arguments{
  \item{x}{SpatVector}
}

\value{
SpatVector
}

\seealso{
\code{\link{sharedPaths}}, \code{\link{topology}}, and 
\code{\link{fillHoles}} to get or remove polygon holes
}


\examples{
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
h <- convHull(v[-12], "NAME_1")
g <- gaps(h)
}

\keyword{methods}
\keyword{spatial}

