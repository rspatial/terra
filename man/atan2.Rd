\name{atan2}

\alias{atan2,SpatRaster,SpatRaster-method}
\alias{atan2}

\title{Two argument arc-tangent}

\description{
For SpatRasters x and y, atan2(y, x) returns the angle in radians for the tangent y/x, handling the case when x is zero. See \code{\link[base]{Trig}}

See \code{\link{Math-methods}} for other trigonometric and mathematical functions that can be used with SpatRasters.
}

\usage{
atan2(y, x) 
}

\arguments{
  \item{y}{SpatRaster}
  \item{x}{SpatRaster}
}

\seealso{ \code{\link{Math-methods}} }


\examples{
r1 <- r2 <- rast(nrow=10, ncol=10)
values(r1) <- (runif(ncell(r1))-0.5) * 10
values(r2) <- (runif(ncell(r1))-0.5) * 10
atan2(r1, r2)
}

\keyword{math}
\keyword{spatial}
