\name{atan2}

\alias{atan2,SpatRaster,SpatRaster-method}
\alias{atan2}
\alias{atan_2,SpatRaster,SpatRaster-method}
\alias{atan_2}

\title{Two argument arc-tangent}

\description{
For SpatRasters x and y, atan2(y, x) returns the angle in radians for the tangent y/x, handling the case when x is zero. See \code{\link[base]{Trig}}

See \code{\link{Math-methods}} for other trigonometric and mathematical functions that can be used with SpatRasters.
}

\usage{
\S4method{atan2}{SpatRaster,SpatRaster}(y, x)

\S4method{atan_2}{SpatRaster,SpatRaster}(y, x, filename, ...)
}

\arguments{
  \item{y}{SpatRaster}
  \item{x}{SpatRaster}
  \item{filename}{character. Output filename}
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}  
}

\seealso{ \code{\link{Math-methods}} }


\examples{
r1 <- rast(nrows=10, ncols=10)
r2 <- rast(nrows=10, ncols=10)
values(r1) <- (runif(ncell(r1))-0.5) * 10
values(r2) <- (runif(ncell(r1))-0.5) * 10
atan2(r1, r2)
}

\keyword{math}
\keyword{spatial}

