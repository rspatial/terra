\name{summary}

\docType{methods}

\alias{summary,SpatRaster-method}

\title{Summary}

\description{
Summarize the values of SpatRaster. A sample is used for very large files. 
}

\usage{
\S4method{summary}{SpatRaster}(object, size=100000, ...)
}

\arguments{
  \item{object}{SpatRaster}
  \item{size}{positive integer. Size of a regular sample used for large datasets}
  \item{...}{additional arguments. None implemented}  
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

