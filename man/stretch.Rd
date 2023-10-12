\name{stretch}

\alias{stretch}
\alias{stretch,SpatRaster-method}
  
\title{Stretch}

\description{
Linear or histogram equalization stretch of values in a SpatRaster. 

For linear stretch, provide the desired output range (\code{minv} and \code{maxv}) and the lower and upper bounds in the original data, either as quantiles (\code{minq} and \code{maxq}, or as cell values (\code{smin} and \code{smax}). If \code{smin} and \code{smax} are both not \code{NA}, \code{minq} and \code{maxq} are ignored.

For histogram equalization, these arguments are ignored, but you can provide the desired scale of the output and the maximum number of cells that is used to compute the histogram (empirical cumulative distribution function). 
}

\usage{
\S4method{stretch}{SpatRaster}(x, minv=0, maxv=255, minq=0, maxq=1, smin=NA, smax=NA,
		histeq=FALSE, scale=1, maxcell=500000, filename="", ...)
}

\arguments{
\item{x}{SpatRaster} 
\item{minv}{numeric >= 0 and smaller than maxv. lower bound of stretched value}
\item{maxv}{numeric <= 255 and larger than maxv. upper bound of stretched value}
\item{minq}{numeric >= 0 and smaller than maxq. lower quantile bound of original value. Ignored if smin is supplied}
\item{maxq}{numeric <= 1 and larger than minq. upper quantile bound of original value. Ignored if smax is supplied}
\item{smin}{numeric < smax. user supplied lower value for the layers, to be used instead of a quantile computed by the function itself}
\item{smax}{numeric > smin. user supplied upper value for the layers, to be used instead of a quantile computed by the function itself}

\item{histeq}{logical. If \code{TRUE} histogram equalization is used instead of linear stretch}
\item{scale}{numeric. The scale (maximum value) of the output if \code{histeq=TRUE}}

\item{maxcell}{positive integer. The size of the regular sample used to compute the histogram}

\item{filename}{character. Output filename}
\item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
}
 
 
\value{
SpatRaster
}
  
\examples{
r <- rast(nc=10, nr=10)
values(r) <- rep(1:25, 4)
rs <- stretch(r)
s <- c(r, r*2)
sr <- stretch(s)
}

\keyword{spatial}

