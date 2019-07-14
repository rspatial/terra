\name{freq}

\docType{methods}

\alias{freq}
\alias{freq,SpatRaster-method}

\title{Frequency table}

\description{
Frequency table of the values of a SpatRaster. \code{NA}s are not counted. 
}

\usage{
\S4method{freq}{SpatRaster}(x, bylayer=TRUE, ...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{bylayer}{logical. If \code{TRUE} tabulation is done by layer}
  \item{...}{additional arguments (none implemented)}
}

\value{
matrix with 2 columns (value, count) or, if \code{bylayer=TRUE} three layers (layer, value, count).
}

\seealso{\code{\link{freq}}}

\examples{
r <- rast(nrow=10, ncol=10)
set.seed(2)
values(r) <- sample(5, ncell(r), replace=TRUE)

freq(r, FALSE)
freq(r)

x <- c(r, r/2)
freq(x, FALSE)
freq(x)

freq(round(x))
}

\keyword{spatial}
\keyword{univar}
