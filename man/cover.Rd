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
\S4method{cover}{SpatRaster,SpatRaster}(x, y, values=NA, filename="", overwrite=FALSE, wopt=list(), ...)

\S4method{cover}{SpatVector,SpatVector}(x, y, identity=FALSE, ...)

}

\arguments{
  \item{x}{SpatRaster}
  \item{y}{SpatRaster}
  \item{values}{numeric. The cell values in \code{x} to be replaced by the values in \code{y}}
  \item{filename}{character. Output filename. Optional}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list. Options for writing files as in \code{\link{writeRaster}}}
  \item{...}{additional arguments. None implemented}  
  \item{identity}{logical. If \code{TRUE} overlapping areas are intersected rather than replaced}
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
cv1 <- cover(p, e)
cv2 <- cover(p, e, identity=TRUE)
}


\keyword{methods}
\keyword{spatial}
