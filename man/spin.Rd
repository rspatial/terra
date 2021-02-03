\name{spin}

\docType{methods}

\alias{spin}
\alias{spin,SpatVector-method}

\title{spin a SpatVector}

\description{
Spin (rotate) the geometry of a SpatVector. 
}

\usage{
\S4method{spin}{SpatVector}(x, angle, x0, y0)
}

\arguments{
 \item{x}{SpatVector}
 \item{angle}{numeric. Angle of rotation in degrees}
 \item{x0}{numeric. x-coordinate of the center of rotation. If missing, the center of the extent of \code{x} is used}
 \item{y0}{numeric. y-coordinate of the center of rotation. If missing, the center of the extent of \code{x} is used}
}

\value{
SpatVector
}

\seealso{\code{\link{rescale}}, \code{\link{t}}, \code{\link{shift}}}

\examples{
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
w <- spin(v, 180)
plot(v)
lines(w, col="red")

# lower-right corner as center
e <- as.vector(ext(v))
x <- spin(v, 45, e[1], e[3])
}

\keyword{spatial}
