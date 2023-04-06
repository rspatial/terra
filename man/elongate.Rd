\name{elongate}

\docType{methods}

\alias{elongate}
\alias{elongate,SpatVector-method}

\title{ elongate lines }

\description{
Elongate SpatVector lines 
}

\usage{
\S4method{elongate}{SpatVector}(x, length=1, flat=FALSE)
}

\arguments{
  \item{x}{SpatVector}
  \item{length}{positive number indicating how much the lines should be elongated at each end. The unit is meter is the crs is lonlat and it is the same as the linear unit of the crs on other cases (also meter in most cases)}
  \item{flat}{logical. If \code{TRUE}, the earth's curvature is ignored for lonlat data, and the distance unit is degrees, not meter}  
}


\value{
SpatVector
}

\seealso{
\code{\link{buffer}}, \code{\link{crop}} and \code{\link{erase}} 
}

\examples{
v <- vect(cbind(c(0,1,2), c(0,0,2)), "lines", crs="lonlat")
e <- elongate(v, 100000)
plot(e)
points(e)
geom(e)
}

\keyword{methods}
\keyword{spatial}

