\name{erase}

\docType{methods}

\alias{erase}
\alias{erase,SpatVector,SpatVector-method}
\alias{erase,SpatVector,SpatExtent-method}

\title{ Erase parts of a SpatVector object}

\description{
Erase parts of a SpatVector with another SpatVector (or SpatExtent). 
The inverse of this can be done with \code{\link{intersect}} and \code{\link{crop}}.
}

\usage{
\S4method{erase}{SpatVector,SpatVector}(x, y)

\S4method{erase}{SpatVector,SpatExtent}(x, y)
}

\arguments{
  \item{x}{SpatVector}
  \item{y}{SpatVector or SpatExtent}
}

\value{
SpatVector or SpatExtent
}

\seealso{The equivalent for SpatRaster is \code{\link{mask}}}


\examples{
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
e <- ext(5.6, 6, 49.55, 49.7)
x <- erase(v, e)

p <- vect("POLYGON ((5.8 49.8, 6 49.9, 6.15 49.8, 6 49.6, 5.8 49.8))")
y <- erase(v, p)
}

\keyword{methods}
\keyword{spatial}

