\name{stretch}

\alias{stretch}
\alias{stretch,SpatRaster-method}
  
\title{Stretch}

\description{
Linear stretch of values in a SpatRaster object. Provide the desired output range (minv and maxv) and the lower and upper bounds in the original data, either as quantiles (\code{minq} and \code{maxq}, or as cell values (\code{smin} and \code{smax}). If \code{smin} and \code{smax} are both not \code{NA}, \code{minq} and \code{maxq} are ignored.
}

\usage{
\S4method{stretch}{SpatRaster}(x, minv=0, maxv=255, minq=0, maxq=1, smin=NA, smax=NA,
          filename="", overwrite=FALSE, wopt=list(), ...)
}

\arguments{
\item{x}{SpatRaster} 
\item{minv}{numeric >= 0 and smaller than maxv. lower bound of stretched value}
\item{maxv}{numeric <= 255 and larger than maxv. upper bound of stretched value}
\item{minq}{numeric >= 0 and smaller than maxq. lower quantile bound of original value. Ignored if smin is supplied}
\item{maxq}{numeric <= 1 and larger than minq. upper quantile bound of original value. Ignored if smax is supplied}
\item{smin}{numeric < smax. user supplied lower value for the layers, to be used instead of a quantile computed by the function itself}
\item{smax}{numeric > smin. user supplied upper value for the layers, to be used instead of a quantile computed by the function itself}

\item{filename}{character. Output filename. Optional}
\item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
\item{wopt}{list. Options for writing files as in \code{\link{writeRaster}}}
\item{...}{ additional arguments as for \code{\link{writeRaster}}}
}
 
 
\value{
SpatRaster
}
  
\examples{
r <- rast(nc=10, nr=10)
values(r) <- rep(1:25, 4)
rs <- stretch(r)
s <- stack(r, r*2)
sr <- stretch(s)
}

\keyword{spatial}

