\name{cartogram}

\docType{methods}

\alias{cartogram}
\alias{cartogram,SpatVector-method}

\title{Cartogram}

\description{ 
Make a cartogram, that is, a map where the area of polygons is made proportional to another variable. This can be a good way to map raw count data (e.g. votes).
}

\usage{
\S4method{cartogram}{SpatVector}(x, var, type="nc")
}

\arguments{
  \item{x}{SpatVector}
  \item{var}{character. A variable name in \code{x}}
  \item{type}{character. Cartogram type, one of "nc" (non-contiguous) or "circles" (dorling)}  
}

\value{
SpatVector
}

\seealso{ \code{\link{plot}}, \code{\link{rescale}} }


\examples{
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
v$value <- 1:12
pnc <- cartogram(v, "value", "nc")
pcirc <- cartogram(v, "value", "circles")
plot(v, col="light gray", border="gray")
lines(pnc, col="red", lwd=2)
lines(pcirc, col="blue", lwd=2)

}


\keyword{methods}
\keyword{spatial}
