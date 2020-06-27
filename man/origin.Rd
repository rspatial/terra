\name{origin}

\alias{origin}
\alias{origin,SpatRaster-method}
  
 
\title{Origin}

\description{
Origin returns the coordinates of the point of origin of a SpatRaster object. This is the point closest to (0, 0) that you could get if you moved towards that point in steps of the x and y resolution.
}

\usage{
origin(x, ...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{...}{additional arguments. None implemented}  
}


\value{
A vector of two numbers (x and y coordinates)
}


\examples{
r <- rast(xmin=-0.5, xmax = 9.5, ncols=10)
origin(r)
}

\keyword{spatial}
