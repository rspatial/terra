\name{summary}

\docType{methods}

\alias{summary}
\alias{summary,SpatRaster-method}

\title{summary}

\description{
Compute summary statistics (min, max, mean, and quartiles) for SpatRaster using base \code{\link[base]{summary}} method. A sample is used for very large files. 
}

\usage{
\S4method{summary}{SpatRaster}(object, size=100000, warn=TRUE, ...)
}

\arguments{
  \item{object}{SpatRaster}
  \item{size}{positive integer. Size of a regular sample used for large datasets (see \code{\link{spatSample}})}
  \item{warn}{logical. If \code{TRUE} a warning is given if a sample is used}  
  \item{...}{additional arguments passed on to the base \code{\link[base]{summary} method}}  
}


\value{matrix with (an estimate of) the median, minimum and maximum values, the first and third quartiles, and the number of cells with \code{NA} values}


\seealso{ \code{\link{global}, \link{quantile}} }

\examples{
set.seed(0)
r <- rast(nrow=10, ncol=10, nlyr=3)
values(r) <- runif(size(r))
summary(r)
}

\keyword{methods}
\keyword{spatial}

