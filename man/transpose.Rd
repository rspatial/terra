\name{transpose}

\docType{methods}

\alias{t}
\alias{t,SpatRaster-method}
\alias{t,SpatVector-method}
\alias{trans}
\alias{trans,SpatRaster-method}
  
\title{Transpose}

\description{
Transpose a SpatRaster or SpatVector
}

\usage{
\S4method{t}{SpatRaster}(x)

\S4method{t}{SpatVector}(x)

\S4method{trans}{SpatRaster}(x, filename="", ...)
}

\arguments{
 \item{x}{SpatRaster or SpatVector}
  \item{filename}{character. Output filename}
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
}

\value{
SpatRaster
}

\seealso{ \code{\link{flip}, \link{rotate}} }



\examples{
r <- rast(nrows=18, ncols=36)
values(r) <- 1:ncell(r)
tr1 <- t(r)
tr2 <- trans(r)
ttr <- trans(tr2)
}

\keyword{spatial}
