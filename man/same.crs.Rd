\name{same.crs}

\alias{same.crs}

  
\title{Compare coordinate reference systems}

\description{
The function takes two coordinate reference system descriptions and compares them for equality.
}

\usage{
same.crs(x, y)
}

\arguments{
\item{x}{character, SpatRaster, SpatVector, CRS, or other object that returns something intelligible with\code{crs(x)}}
\item{y}{same types as for \code{x}}
}

\value{
logical
}

\examples{
r <- rast()
same.crs(r, "+proj=longlat")

same.crs(r, "+proj=utm +zone=1")
}

\keyword{spatial}

