\name{flip}

\docType{methods}

\alias{flip}
\alias{flip,SpatRaster-method}
  
\title{Flip a raster}

\description{
Flip the values of a SpatRaster by inverting the order of the rows (\code{vertical=TRUE}) or the columns (\code{vertical=FALSE}).
}

\usage{
\S4method{flip}{SpatRaster}(x, vertical=TRUE, filename="", overwrite=FALSE, wopt=list(), ...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{vertical}{logical. If \code{TRUE}, values are flipped by rows, if \code{FALSE}, values are flipped by columns}
  \item{filename}{character. Output filename. Optional}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list. Options for writing files as in \code{\link{writeRaster}}}
  \item{...}{additional arguments. None implemented}  
}

\value{
SpatRaster
}


\seealso{ \code{\link{transpose}}, \code{\link{rotate}} }



\examples{
r <- rast(nrow=18, ncol=36)
m <- matrix(1:ncell(r), nrow=18)
values(r) <- as.vector(t(m))
rx <- flip(r, vertical=FALSE)
values(r) <- as.vector(m)
ry <- flip(r, vertical=TRUE)
}

\keyword{spatial}
