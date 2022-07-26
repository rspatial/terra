\name{cover}

\docType{methods}

\alias{cover}
\alias{cover,SpatRaster,SpatRaster-method}
\alias{cover,SpatVector,SpatVector-method}


\title{Replace values with values from another object}

\description{
Replace \code{NA} or other values in SpatRaster \code{x} with the values of SpatRaster \code{y}

For polygons: areas of \code{x} that overlap with \code{y} are replaced by \code{y} or, if \code{identity=TRUE} intersected with \code{y}.

}

\usage{
\S4method{cover}{SpatRaster,SpatRaster}(x, y, values=NA, filename="", ...)

\S4method{cover}{SpatVector,SpatVector}(x, y, identity=FALSE, expand=TRUE)
}

\arguments{
  \item{x}{SpatRaster or SpatVector}
  \item{y}{Same as \code{x}}
  \item{values}{numeric. The cell values in \code{x} to be replaced by the values in \code{y}}
  \item{filename}{character. Output filename}
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
  \item{identity}{logical. If \code{TRUE} overlapping areas are intersected rather than replaced}
  \item{expand}{logical. Should parts of \code{y} that are outside of \code{x} be included?}
}


\value{
SpatRaster
}

\examples{
r1 <- r2 <- rast(ncols=36, nrows=18)
values(r1) <- 1:ncell(r1)
values(r2) <- runif(ncell(r2))
r2 <- classify(r2, cbind(-Inf, 0.5, NA))
r3 <- cover(r2, r1)


p <- vect(system.file("ex/lux.shp", package="terra"))
e <- as.polygons(ext(6, 6.4, 49.75, 50))
values(e) <- data.frame(y=10)

cv <- cover(p, e)
plot(cv, col=rainbow(12))
ci <- cover(p, e, identity=TRUE)
lines(e, lwd=3)

plot(ci, col=rainbow(12))
lines(e, lwd=3)

}


\keyword{methods}
\keyword{spatial}
