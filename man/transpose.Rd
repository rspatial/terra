\name{transpose}

\docType{methods}

\alias{t}
\alias{t,SpatRaster-method}
\alias{transpose}
\alias{transpose,SpatRaster-method}
  
\title{Transpose}

\description{
Transpose a SpatRaster
}

\usage{
\S4method{t}{SpatRaster}(x)

\S4method{transpose}{SpatRaster}(x, filename="", overwrite=FALSE, wopt=list(), ...)
}

\arguments{
 \item{x}{SpatRaster}
  \item{filename}{character. Output filename. Optional}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list. Options for writing files as in \code{\link{writeRaster}}}
  \item{...}{additional arguments. None implemented}  
 
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
