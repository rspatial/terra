\name{freq}

\docType{methods}

\alias{freq}
\alias{freq,SpatRaster-method}

\title{Frequency table}

\description{
Frequency table of the values of a SpatRaster. \code{NA}s are not counted. 
}

\usage{
\S4method{freq}{SpatRaster}(x, digits=0, value=NULL, bylayer=TRUE)
}

\arguments{
  \item{x}{SpatRaster}
  \item{bylayer}{logical. If \code{TRUE} tabulation is done by layer}
  \item{digits}{integer. Used for rounding the values before tabulation. Ignored if \code{NA}}
  \item{value}{numeric. An optional single value to only count the number of cells with that value}
}

\value{
matrix with 2 columns (value, count) or, if \code{bylayer=TRUE} three layers (layer, value, count).
}

\seealso{\code{\link{freq}}}

\examples{
r <- rast(nrows=10, ncols=10)
set.seed(2)
values(r) <- sample(5, ncell(r), replace=TRUE)

freq(r, FALSE)
freq(r)

x <- c(r, r/3)
freq(x, FALSE)
freq(x)
freq(x, digits=1)
freq(x, digits=-1)

freq(x, value=5)

}

\keyword{spatial}
\keyword{univar}
