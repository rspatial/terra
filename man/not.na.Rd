\name{not.na}

\docType{methods}

\alias{not.na}
\alias{not.na,SpatRaster-method}


\title{is not NA}

\description{
Shortcut method to avoid the two-step \code{!is.na(x)}
}

\usage{
\S4method{not.na}{SpatRaster}(x, falseNA, filename="", ...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{falseNA}{logical. Should the result be \code{TRUE, NA} instead of \code{TRUE, FALSE}?}
  \item{filename}{character. Output filename}
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
}


\seealso{
\code{\link{Compare-methods}}
}

\value{
SpatRaster
}

\examples{
r <- rast(ncols=10, nrows=10, vals=1)
r[10:20] <- NA
x <- not.na(r)
}

\keyword{methods}
\keyword{spatial}

