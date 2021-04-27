\name{flip}

\docType{methods}

\alias{flip}
\alias{flip,SpatRaster-method}
\alias{flip,SpatVector-method}

\alias{rev}
\alias{rev,SpatRaster-method}
  
\title{Flip or reverse a raster}

\description{
Flip the values of a SpatRaster by inverting the order of the rows (\code{vertical=TRUE}) or the columns (\code{vertical=FALSE}).

\code{rev} is the same as a horizontal *and* a vertical flip.
}

\usage{
\S4method{flip}{SpatRaster}(x, direction="vertical", filename="", ...)

\S4method{flip}{SpatVector}(x, direction="vertical")

\S4method{rev}{SpatRaster}(x)
}

\arguments{
  \item{x}{SpatRaster or SpatVector}
  \item{direction}{character. Should (partially) match "vertical" to flip by rows, or "horizontal" to flip by columns}
  \item{filename}{character. Output filename}
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
}

\value{
SpatRaster
}

\seealso{ \code{\link{trans}}, \code{\link{rotate}} }


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
