\name{subst}

\docType{methods}

\alias{subst}
\alias{subst,SpatRaster-method}

\title{replace cell values}

\description{
Substitute(replace) cell values of a SpatRaster with a new value. See \code{\link{classify}} for more complex/flexible replacement. 
}

\usage{
\S4method{subst}{SpatRaster}(x, from, to, filename="", ...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{from}{numeric value(s)}
  \item{to}{a single numeric value} 
  \item{filename}{character. Output filename}
  \item{...}{Additional arguments for writing files as in \code{\link{writeRaster}}}
}

\value{
SpatRaster
}

\seealso{\code{\link{classify}}}

\examples{
r <- rast(ncol=5, nrow=5, xmin=0, xmax=1, ymin=0, ymax=1, crs="")
r <- init(r, 1:6)
x <- subst(r, 3, 7)
x <- subst(r, 2:3, NA)
x <- subst(x, NA, 10)
}

\keyword{spatial}

