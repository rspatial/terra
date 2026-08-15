\name{has.z}

\docType{methods}

\alias{has.z}
\alias{has.z,SpatVector-method}


\title{Does a SpatVector have Z coordinates?}

\description{
Returns \code{TRUE} if any geometry part (or polygon hole) stores Z values.

Z storage is optional: empty Z means a conventional 2D geometry. When Z is
present it is included in \code{\link{geom}} and \code{\link{crds}}, preserved
by \code{\link{writeVector}} / \code{\link{vect}} file I/O (25D geometries),
and by coercion to/from \pkg{sf} (see also
\url{https://github.com/rspatial/terra/issues/824}).
}

\usage{
\S4method{has.z}{SpatVector}(x)
}

\arguments{
  \item{x}{SpatVector}
}

\value{
logical
}

\seealso{
 \code{\link{geom}, \link{crds}, \link{vect}}
}

\examples{
v <- vect(cbind(1:3, 1:3))
has.z(v)
}

\keyword{methods}
\keyword{spatial}
