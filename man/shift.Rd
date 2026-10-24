\name{shift}

\docType{methods}

\alias{shift}
\alias{shift,SpatRaster-method}
\alias{shift,SpatExtent-method}
\alias{shift,SpatVector-method}
  
\title{Shift}

\description{
Shift a SpatRaster, SpatVector or SpatExtent to another location.
}

\usage{
\S4method{shift}{SpatRaster}(x, dx=0, dy=0, filename="", to, ...)

\S4method{shift}{SpatVector}(x, dx=0, dy=0, to)

\S4method{shift}{SpatExtent}(x, dx=0, dy=0, to)
}

\arguments{
  \item{x}{SpatRaster, SpatVector or SpatExtent}
  \item{dx}{numeric. The shift in horizontal direction}
  \item{dy}{numeric. The shift in vertical direction}
  \item{to}{SpatExtent, or an object from which an extent can be extracted (see \code{\link{ext}}); or two coordinates (x, y) for the new center of \code{x}. One of the two coordinates may be \code{NA} to shift only in the other direction. If \code{to} is used, \code{dx} and \code{dy} cannot be used}
  \item{filename}{character. Output filename}
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
}
 
\value{
Same as \code{x}
}


\seealso{ \code{\link{flip}}, \code{\link{rotate}}}

\examples{
r <- rast(xmin=0, xmax=1, ymin=0, ymax=1)
r <- shift(r, dx=1, dy=-1)

e <- ext(r)
shift(e, 5, 5)

# center on another object or on coordinates
e2 <- ext(10, 12, 4, 8)
shift(e, to=e2)
shift(e, to=c(11, 6))
shift(e, to=c(11, NA))
}

\keyword{spatial}
