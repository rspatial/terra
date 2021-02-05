\name{rescale}

\docType{methods}

\alias{rescale}
\alias{rescale,SpatRaster-method}
\alias{rescale,SpatVector-method}

\title{rescale}

\description{
Rescale a SpatVector or SpatRaster. This may be useful to make small \code{\link{inset}} maps. 
}

\usage{
\S4method{rescale}{SpatRaster}(x, f=0.5, x0, y0, ...)

\S4method{rescale}{SpatVector}(x, f=0.5, x0, y0, ...)
}

\arguments{
  \item{x}{SpatVector or SpatRaster}
  \item{f}{numeric. The scaling factor (a fraction)}
  \item{x0}{numeric. x-coordinate of the center of rescaling. If missing, the center of the extent of \code{x} is used}
  \item{y0}{numeric. y-coordinate of the center of rescaling. If missing, the center of the extent of \code{x} is used}
  \item{...}{additional arguments. None implemented}  
}
 
\value{
Same as \code{x}
}


\seealso{\code{\link{t}}, \code{\link{shift}}, \code{\link{flip}}, \code{\link{rotate}}, \code{\link{inset}}}

\examples{
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
w <- rescale(v, 0.2)
plot(v)
lines(w, col="red")
}

\keyword{spatial}
