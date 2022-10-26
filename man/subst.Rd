\name{subst}

\docType{methods}

\alias{subst}
\alias{subst,SpatRaster-method}

\title{replace cell values}

\description{
Substitute(replace) cell values of a SpatRaster with a new value. See \code{\link{classify}} for more complex/flexible replacement. 
}

\usage{
\S4method{subst}{SpatRaster}(x, from, to, others=NULL, raw=FALSE, filename="", ...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{from}{numeric value(s). Normally a vector of the same length as `to`. If \code{x} has multiple layers, it can also be a matrix of numeric value(s) where \code{nrow(x) == length(to)}. In that case the output has a single layer, with values based on the combination of the values of the input layers}
  \item{to}{numeric value(s). Normally a vector of the same length as `from`. If \code{x} has a single layer, it can also be a matrix of numeric value(s) where \code{nrow(x) == length(from)}. In that case the output has multiple layers, one for each column in \code{to}} 
  \item{others}{numeric. If not \code{NULL} all values that are not matched are set to this value. Otherwise they retain their original value.} 
  \item{raw}{logical. If \code{TRUE}, the values in from and to are the raw cell values, not the categorical labels. Only relevant if \code{is.factor(x)}}
  \item{filename}{character. Output filename}
  \item{...}{Additional arguments for writing files as in \code{\link{writeRaster}}}
}

\value{
SpatRaster
}

\seealso{\code{\link{classify}}}

\examples{
r <- rast(ncols=5, nrows=5, xmin=0, xmax=1, ymin=0, ymax=1, crs="")
r <- init(r, 1:6)
x <- subst(r, 3, 7)
x <- subst(r, 2:3, NA)
x <- subst(x, NA, 10)

# multiple output layers
z <- subst(r, 2:3, cbind(20,30))

# multiple input layers
rr <- c(r, r+1, r+2)
m <- rbind(c(1:3), c(3:5))
zz <- subst(rr, m, c(100, 200))
}

\keyword{spatial}

