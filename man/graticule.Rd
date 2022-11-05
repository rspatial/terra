\name{graticule}

\docType{methods}

\alias{graticule}

\title{ 
Create a graticule
}

\description{
Create a graticule. That is, a grid of lon/lat lines that can be used to on a projected map.

The object returned, a SpatGraticule, can be plotted with \code{plot} and \code{lines}. There is also a \code{crop} method.
}

\usage{
graticule(lon=30, lat=30, crs="+proj=longlat")
}

\arguments{
  \item{lon}{numeric. Either single number (the interval between longitudes) or a vector with longitudes}
  \item{lat}{numeric. Either single number (the interval between latitudes) or a vector with latitudes}
  \item{crs}{character}
} 

\value{
SpatGraticule
}


\seealso{
\code{\link[terra:plot_graticule]{plot<SpatGraticule>}}.
}

\examples{
g <- graticule(60, 30, crs="+proj=robin")
g

graticule(90, c(-90, -60, -23.5, 0, 23.5, 60, 90), crs="+proj=robin")

}

\keyword{methods}
\keyword{spatial}

