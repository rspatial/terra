\name{normalize.longitude}

\docType{methods}

\alias{normalize.longitude}
\alias{normalize.longitude,SpatVector-method}

\title{normalize vector data that crosses the dateline}

\description{
Normalize the longitude of geometries, move them if they are outside of the -180 to 180 degrees range. 
} 

\usage{
\S4method{normalize.longitude}{SpatVector}(x)
}

\arguments{
  \item{x}{SpatVector}
}

\value{SpatVector}

\seealso{\link{rotate} for SpatRaster}


\examples{
p <- vect("POLYGON ((120 10, 230 75, 230 -75, 120 10))")
normalize.longitude(p)
}

\keyword{methods}
\keyword{spatial}

