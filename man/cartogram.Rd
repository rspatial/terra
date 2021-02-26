\name{cartogram}

\docType{methods}

\alias{cartogram}
\alias{cartogram,SpatVector-method}

\title{Cartogram}

\description{ 
Make a cartogram, that is, a map where the area of polygons is made proportional to another variable. This can be a good way to map raw count data (e.g. votes).
}

\usage{
\S4method{cartogram}{SpatVector}(x, var, type)
}

\arguments{
  \item{x}{SpatVector}
  \item{var}{character. A variable name in \code{x}}
  \item{type}{character. Cartogram type, only "nc" (non-contiguous) is currently supported}  
}

\value{
SpatVector
}

\seealso{ \code{\link{plot}}, \code{\link{rescale}} }


\examples{
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
v$value <- 1:12
p <- cartogram(v, "value", "nc")
plot(v, col="light gray", border="gray")
lines(p, col="red", lwd=2)
}


\keyword{methods}
\keyword{spatial}
