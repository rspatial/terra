\name{summarize-generics}

\docType{methods}

\alias{Summary-methods}

\alias{median}
\alias{mean}
\alias{min}
\alias{max}
\alias{range}
\alias{which.min}
\alias{which.max}

\alias{Summary,SpatExtent-method}
\alias{Summary,SpatRaster-method}
\alias{mean,SpatRaster-method}
\alias{mean,SpatExtent-method}
\alias{median,SpatRaster-method}
\alias{stdev}
\alias{stdev,SpatRaster-method}
\alias{min,SpatRaster-method}
\alias{max,SpatRaster-method}
\alias{range,SpatRaster-method}
\alias{which.min,SpatRaster-method}
\alias{which.max,SpatRaster-method}


\title{Summarize}

\description{

Compute summary statistics for cells, either across layers or between layers (parallel summary). 

The following summary methods are available for SpatRaster: \code{any, all, max, min, mean, median, prod, range, stdev, sum, which.min, which.max}. See \code{\link{modal}} to compute the mode and \code{\link{app}} to compute summary statistics that are not included here.
 
Because generic functions are used, the method applied is chosen based on the first argument: "\code{x}". This means that if \code{r} is a SpatRaster, \code{mean(r, 5)} will work, but \code{mean(5, r)} will not work.

The \code{mean} method has an argument "trim" that is ignored.

The \code{stdev} method returns the population standard deviation, computed as:

\code{f <- function(x) sqrt(sum((x-mean(x))^2) / length(x))}

This is different than the sample standard deviation returned by \code{sd} (which uses \code{n-1} as denominator). Function \code{f} above is equivalent to function \code{g} below

\code{g <- function(x) sqrt(sum((x-mean(x))^2) / length(x))}
}


\usage{
\S4method{min}{SpatRaster}(x, ..., na.rm=FALSE)

\S4method{max}{SpatRaster}(x, ..., na.rm=FALSE)

\S4method{range}{SpatRaster}(x, ..., na.rm=FALSE)

\S4method{mean}{SpatRaster}(x, ..., trim=NA, na.rm=FALSE)

\S4method{median}{SpatRaster}(x, na.rm=FALSE)

\S4method{stdev}{SpatRaster}(x, ..., na.rm=FALSE)

\S4method{which.min}{SpatRaster}(x)

\S4method{which.max}{SpatRaster}(x)
}


\arguments{
  \item{x}{SpatRaster}
  \item{...}{additional SpatRaster objects or numeric values}
  \item{trim}{ignored}
  \item{na.rm}{logical. If \code{TRUE}, \code{NA} values are ignored. If \code{FALSE}, \code{NA} is returned if \code{x} has any \code{NA} values}
}


\value{SpatRaster}

\seealso{\code{\link{app}}, \code{\link{Math-methods}}, \code{\link{modal}}, \code{\link{which.lyr}} }


\examples{
set.seed(0)
r <- rast(nrow=10, ncol=10, nlyr=3)
values(r) <- runif(size(r))

x <- mean(r)
# note how this returns one layer
x <- sum(c(r, r[[2]]), 5)

# and this returns three layers
y <- sum(r, r[[2]], 5)

max(r)
max(r, 0.5)

y <- stdev(r)
# not the same as 
yy <- app(r, sd)

z <- stdev(r, r*2)
}

\keyword{methods}
\keyword{spatial}

