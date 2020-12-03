\name{rescale}

\docType{methods}

\alias{rescale}
\alias{rescale,SpatVector-method}
  
\title{rescale}

\description{
rescale a SpatVector
}

\usage{
\S4method{rescale}{SpatVector}(x, f, ...)
}

\arguments{
  \item{x}{SpatVector}
  \item{f}{numeric. The scaling factor (a fraction)}
  \item{...}{additional arguments. None implemented}  
}
 
\value{
SpatVector
}


\seealso{ \code{\link{t}}, \code{\link{shift}}, \code{\link{flip}}, \code{\link{rotate}}}

\examples{
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
w <- rescale(v, 0.2)
plot(v)
lines(w, col="red")
}

\keyword{spatial}
