\name{origin}

\alias{origin}
\alias{origin,SpatRaster-method}
\alias{origin<-}
\alias{origin<-,SpatRaster-method}

 
\title{Origin}

\description{
Get or set the coordinates of the point of origin of a SpatRaster. This is the point closest to (0, 0) that you could get if you moved towards that point in steps of the x and y resolution.
}

\usage{
\S4method{origin}{SpatRaster}(x)

\S4method{origin}{SpatRaster}(x)<-value
}

\arguments{
  \item{x}{SpatRaster}
  \item{value}{numeric vector of length 1 or 2}
}


\value{
A vector of two numbers (x and y coordinates)
}


\examples{
r <- rast(xmin=-0.5, xmax = 9.5, ncols=10)
origin(r)
origin(r) <- c(0,0)
r
}

\keyword{spatial}
