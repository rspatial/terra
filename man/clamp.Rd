\name{clamp}

\alias{clamp}
\alias{clamp,SpatRaster-method}
\alias{clamp,numeric-method}

\title{Clamp values}

\description{
Clamp values to a minimum and maximum value. That is, all values below a lower threshold value and above the upper threshold value become either \code{NA}, or, if \code{values=TRUE}, become the threshold value  
}

\usage{
\S4method{clamp}{SpatRaster}(x, lower=-Inf, upper=Inf, values=TRUE, filename="", ...)

\S4method{clamp}{numeric}(x, lower=-Inf, upper=Inf, values=TRUE, ...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{lower}{numeric with the lowest acceptable value (you can specify a different value for each layer). Or a SpatRaster that has a single layer or the same number of layers as \code{x}}
  \item{upper}{numeric with the highest acceptable value (you can specify a different value for each layer). Or a SpatRaster that has a single layer or the same number of layers as \code{x}}
  \item{values}{logical. If \code{FALSE} values outside the clamping range become \code{NA}, if \code{TRUE}, they get the extreme values}
  \item{filename}{character. Output filename}
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
}

\value{
SpatRaster
}

\seealso{
\code{\link{classify}}
}

\examples{
r <- rast(ncols=10, nrows=10)
values(r) <- 1:ncell(r)
rc <- clamp(r, 25, 75) 
rc
}

\keyword{spatial}

