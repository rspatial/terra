\name{local}

\docType{methods}

\alias{Summary-methods}

\alias{Summary,SpatRaster-method}
\alias{mean,SpatRaster-method}
\alias{stdev}
\alias{stdev,SpatRaster-method}
\alias{min,SpatRaster-method}
\alias{max,SpatRaster-method}
\alias{range,SpatRaster-method}

\title{Local statistics}

\description{

Compute cell (pixel) level "local" statistics across layers, and most return a single layer (the exception is \code{range}). 

The following summary methods are available for SpatRaster: \code{any, all, max, min, mean, prod, range, stdev, sum}

See \code{app} to compute statistics that are not included here (e.g. median) or any other custom functions.
 
See \code{\link{modal}} to compute the mode.  
 
Because generic functions are used, the method applied is chosen based on the first argument: "\code{x}". This means that if \code{r} is a SpatRaster, \code{mean(r, 5)} will work, but \code{mean(5, r)} will not work.

The \code{mean} method has an argument "trim" that is ignored.

The \code{stdev} mehtod returns the population standard deviation, computed as such:

\code{f <- function(x) sqrt(sum((x-mean(x))^2) / length(x))}

This is different than the sample standard deviation retunred by \code{sd} (which uses \code{n-1} as denominator). Function \code{f} above is equivalent to function \code{g} below

\code{g <- function(x) sqrt(sum((x-mean(x))^2) / length(x))}
}


\usage{
\S4method{min}{SpatRaster}(x, ..., na.rm=FALSE)
\S4method{max}{SpatRaster}(x, ..., na.rm=FALSE)
\S4method{range}{SpatRaster}(x, ..., na.rm=FALSE)
\S4method{mean}{SpatRaster}(x, ..., trim=NA, na.rm=FALSE)
\S4method{stdev}{SpatRaster}(x, ..., na.rm=FALSE)
}


\arguments{
  \item{x}{SpatRaster}
  \item{...}{additional argument of the same type as \code{x} or numeric}
  \item{trim}{ignored}
  \item{na.rm}{logical. If \code{TRUE}, \code{NA} values are ignored. If \code{FALSE}, \code{NA} is returned if \code{x} has any \code{NA} values}
}


\value{SpatRaster}

\seealso{ \code{\link{Math-methods}, \link{modal}} }


\examples{
set.seed(0)
r <- rast(nrow=10, ncol=10, nlyr=3)
values(r) <- runif(size(r))

x <- mean(r)
x <- sum(r, r[[2]], 3)
max(r)
max(r, 0.5)

y <- stdev(r)
# not the same as 
yy <- app(r, sd)

z <- stdev(r, r*2)
}

\keyword{methods}
\keyword{spatial}

