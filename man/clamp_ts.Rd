\name{clamp_ts}

\alias{clamp_ts}
\alias{clamp_ts,SpatRaster-method}
\alias{clamp_ts,numeric-method}

\title{clamp time series data}

\description{
clamp time-series datat that are S shaped. The value in layers before the minimum value in a cell can be set to that minimum value, and the value in layers after the maximum value for a cell can be set to that maximum value.   
}

\usage{
\S4method{clamp_ts}{SpatRaster}(x, min=FALSE, max=TRUE, filename="", ...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{min}{logical. If \code{TRUE} the time-series is clamped to the minimum value}
  \item{max}{logical. If \code{TRUE} the time-series is clamped to the maximum value}
  \item{filename}{character. Output filename}
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
}

\value{
SpatRaster
}

\seealso{
\code{\link{clamp}}, \code{\link{cummin}}, \code{\link{cummax}}
}

\examples{
sigm <- function(x) { .8 / (1 + exp(-(x-10))) + runif(length(x))/4 }
r <- rast(ncols=10, nrows=10, nlyr=50)
s <- seq(5.2, 15,.2)
set.seed(1)
values(r) <- t(replicate(100, sigm(s)))

x <- clamp_ts(r, TRUE, TRUE) 

plot(unlist(r[4]))
lines(unlist(x[4]))

}

\keyword{spatial}

