\name{transpose}

\docType{methods}

\alias{t}
\alias{t,SpatRaster-method}
\alias{t,SpatVector-method}
\alias{transpose}
\alias{transpose,SpatRaster-method}
  
\title{Transpose}

\description{
Transpose a SpatRaster
}

\usage{
\S4method{t}{SpatRaster}(x)

\S4method{t}{SpatVector}(x)

\S4method{transpose}{SpatRaster}(x, filename="", ...)
}

\arguments{
 \item{x}{SpatRaster}
  \item{filename}{character. Output filename}
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
}

\value{
SpatRaster
}

\seealso{ \code{\link{flip}, \link{rotate}} }



\examples{
r <- rast(nrow=18, ncol=36)
values(r) <- 1:ncell(r)
tr1 <- t(r)
tr2 <- transpose(r)
ttr <- transpose(tr2)
}

\keyword{spatial}
