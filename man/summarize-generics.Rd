\name{summarize}

\docType{methods}

\alias{Summary-methods}

\alias{median}
\alias{mean}
\alias{min}
\alias{max}
\alias{which.min}
\alias{which.max}
\alias{any}
\alias{all}
\alias{prod}
\alias{range}
\alias{sum}
\alias{anyNA}
\alias{allNA}
\alias{noNA}
\alias{countNA}

\alias{mean,SpatRaster-method}
\alias{mean,SpatVector-method}
\alias{mean,SpatExtent-method}
\alias{median,SpatRaster-method}
\alias{median,SpatVector-method}
\alias{min,SpatRaster-method}
\alias{max,SpatRaster-method}
\alias{range,SpatRaster-method}
\alias{which.min,SpatRaster-method}
\alias{which.max,SpatRaster-method}

\alias{any,SpatRaster-method}
\alias{all,SpatRaster-method}
\alias{prod,SpatRaster-method}
\alias{range,SpatRaster-method}
\alias{sum,SpatRaster-method}
\alias{anyNA,SpatRaster-method}
\alias{allNA,SpatRaster-method}
\alias{noNA,SpatRaster-method}
\alias{countNA,SpatRaster-method}

\alias{stdev}
\alias{stdev,SpatRaster-method}


\title{Summarize}

\description{

Compute summary statistics for cells, either across layers or between layers (parallel summary). 

The following summary methods are available for SpatRaster: \code{any, anyNA, all, allNA, max, min, mean, median, prod, range, stdev, sum, which.min, which.max}. See \code{\link{modal}} to compute the mode and \code{\link{app}} to compute summary statistics that are not included here.
 
Because generic functions are used, the method applied is chosen based on the first argument: "\code{x}". This means that if \code{r} is a SpatRaster, \code{mean(r, 5)} will work, but \code{mean(5, r)} will not work.

The \code{mean} method has an argument "trim" that is ignored.

If \code{pop=TRUE} \code{stdev} computes the population standard deviation, computed as:

\code{f <- function(x) sqrt(sum((x-mean(x))^2) / length(x))}

This is different than the sample standard deviation returned by \code{sd} (which uses \code{n-1} as denominator). 
}


\usage{
\S4method{min}{SpatRaster}(x, ..., na.rm=FALSE)

\S4method{max}{SpatRaster}(x, ..., na.rm=FALSE)

\S4method{range}{SpatRaster}(x, ..., na.rm=FALSE)

\S4method{prod}{SpatRaster}(x, ..., na.rm=FALSE)

\S4method{sum}{SpatRaster}(x, ..., na.rm=FALSE)

\S4method{any}{SpatRaster}(x, ..., na.rm=FALSE)

\S4method{all}{SpatRaster}(x, ..., na.rm=FALSE)

\S4method{range}{SpatRaster}(x, ..., na.rm=FALSE)

\S4method{which.min}{SpatRaster}(x)

\S4method{which.max}{SpatRaster}(x)

\S4method{stdev}{SpatRaster}(x, ..., pop=TRUE, na.rm=FALSE)

\S4method{mean}{SpatRaster}(x, ..., trim=NA, na.rm=FALSE)

\S4method{median}{SpatRaster}(x, na.rm=FALSE, ...)

\S4method{anyNA}{SpatRaster}(x)

\S4method{countNA}{SpatRaster}(x, n=0)

\S4method{noNA}{SpatRaster}(x, falseNA=FALSE)

\S4method{allNA}{SpatRaster}(x, falseNA=FALSE)
}


\arguments{
  \item{x}{SpatRaster}
  \item{...}{additional SpatRaster objects or numeric values; and arguments \code{par} for parallel summarization (see Details), and \code{filename}, \code{overwrite} and \code{wopt} as for \code{\link{writeRaster}}}
  \item{na.rm}{logical. If \code{TRUE}, \code{NA} values are ignored. If \code{FALSE}, \code{NA} is returned if \code{x} has any \code{NA} values}
  \item{trim}{ignored}
  \item{pop}{logical. If \code{TRUE}, the population standard deviation is computed. Otherwise the sample standard deviation is computed}
  \item{falseNA}{logical. If \code{TRUE}, cells that would otherwise be \code{FALSE} are set to \code{NA}}
  \item{n}{integer. If \code{n > 0}, cell values are \code{TRUE} if at least \code{n} of its layers are \code{NA}}
}


\value{SpatRaster}

\details{
Additional argument \code{par} can be used for "parallel" summarizing a SpatRaster and a numeric or logical value. If a SpatRaster \code{x} has three layers, \code{max(x, 5)} will return a single layer (the number five is treated as a layer in which all cells have value five). In contrast \code{max(x, 5, par=TRUE)} returns three layers (the number five is treated as another SpatRaster with a single layer with all cells having the value five.
}


\seealso{\code{\link{app}}, \code{\link{Math-methods}}, \code{\link{modal}}, \code{\link{which.lyr}} }


\examples{
set.seed(0)
r <- rast(nrows=10, ncols=10, nlyrs=3)
values(r) <- runif(ncell(r) * nlyr(r))

x <- mean(r)
# note how this returns one layer
x <- sum(c(r, r[[2]]), 5)

# and this returns three layers
y <- sum(r, r[[2]], 5)

max(r)

## when adding a number, do you want 1 layer or all layers?
# 1 layer
max(r, 0.5)

# all layers
max(r, 0.5, par=TRUE)

y <- stdev(r)
# not the same as 
yy <- app(r, sd)

z <- stdev(r, r*2)

x <- mean(r, filename=paste0(tempfile(), ".tif"))


v <- values(r)
set.seed(3)
v[sample(length(v), 50)] <- NA
values(r) <- v
is.na(r)
anyNA(r)
allNA(r)
countNA(r)
countNA(r, 2)
}

\keyword{methods}
\keyword{spatial}

