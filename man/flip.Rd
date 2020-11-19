\name{flip}

\docType{methods}

\alias{flip}
\alias{flip,SpatRaster-method}

\alias{rev}
\alias{rev,SpatRaster-method}
  
\title{Flip or reverse a raster}

\description{
Flip the values of a SpatRaster by inverting the order of the rows (\code{vertical=TRUE}) or the columns (\code{vertical=FALSE}).

\code{rev} is the same as a horizontal *and* a vertical flip.
}

\usage{
\S4method{flip}{SpatRaster}(x, vertical=TRUE, filename="", overwrite=FALSE, wopt=list(), ...)

\S4method{rev}{SpatRaster}(x)
}

\arguments{
  \item{x}{SpatRaster}
  \item{direction}{character. Should (partially) match "vertical" to flip by rows, or "horizontal" to flip by columns}
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
rx <- flip(r, direction="h")

values(r) <- as.vector(m)
ry <- flip(r, direction="v")

v <- rev(r)
}

\keyword{spatial}
