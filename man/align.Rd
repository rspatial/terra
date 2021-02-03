\name{align}

\alias{align}
\alias{align,SpatExtent,SpatRaster-method} 
\alias{align,SpatExtent,numeric-method} 
 
 
\title{Align a SpatExtent}

\description{
Align an SpatExtent with a SpatRaster This can be useful to create a new SpatRaster with the same origin and resolution as an existing SpatRaster. Do not use this to force data to match that really does not match (use e.g. \code{\link{resample}} or (dis)aggregate for this).

It is also possible to align a SpatExtent to a clean divisor.
}


\usage{
\S4method{align}{SpatExtent,SpatRaster}(x, y, snap="near")

\S4method{align}{SpatExtent,numeric}(x, y)
}

\arguments{
  \item{x}{SpatExtent}  
  \item{y}{SpatRaster or numeric}
  \item{snap}{Character. One of "near", "in", or "out", to determine in which direction the extent should be aligned. To the nearest border, inwards or outwards}  
}
 
\value{
SpatExtent
}
 
\seealso{ \code{\link{ext}}, \code{\link{draw}} }

\examples{
r <- rast()
e <- ext(-10.1, 9.9, -20.1, 19.9)
ea <- align(e, r)
e
ext(r)
ea

align(e, 0.5)
}

\keyword{spatial}
