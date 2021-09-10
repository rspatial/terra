\name{summary}

\docType{methods}

\alias{summary}
\alias{summary,SpatRaster-method}
\alias{summary,SpatVector-method}
\alias{Summary,SpatExtent-method}
\alias{Summary,SpatRaster-method}
\alias{Summary,SpatVector-method}

\title{summary}

\description{
Compute summary statistics (min, max, mean, and quartiles) for SpatRaster using base \code{\link[base]{summary}} method. A sample is used for very large files. 

For single or other statistics see \code{\link{Summary-methods}}, \code{\link{global}}, and \code{\link{quantile}} 
}

\usage{
\S4method{summary}{SpatRaster}(object, size=100000, warn=TRUE, ...)

\S4method{summary}{SpatVector}(object, ...)
}

\arguments{
  \item{object}{SpatRaster or SpatVector}
  \item{size}{positive integer. Size of a regular sample used for large datasets (see \code{\link{spatSample}})}
  \item{warn}{logical. If \code{TRUE} a warning is given if a sample is used}  
  \item{...}{additional arguments passed on to the base \code{\link[base]{summary} method}}  
}


\value{matrix with (an estimate of) the median, minimum and maximum values, the first and third quartiles, and the number of cells with \code{NA} values}


\seealso{ \code{\link{Summary-methods}}, \code{\link{global}}, \code{\link{quantile}} }

\examples{
set.seed(0)
r <- rast(nrows=10, ncols=10, nlyrs=3)
values(r) <- runif(nlyr(r)*ncell(r))
summary(r)
}

\keyword{methods}
\keyword{spatial}

