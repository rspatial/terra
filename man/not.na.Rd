\name{not.na}

\docType{methods}

\alias{not.na}
\alias{not.na,SpatRaster-method}


\title{is not NA}

\description{
Shortcut method to avoid the two-step \code{!is.na(x)}
}

\usage{
\S4method{not.na}{SpatRaster}(x, falseNA=FALSE, filename="", ...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{falseNA}{logical. If \code{TRUE}, the output values are \code{TRUE} for cells that are not \{NA}, and \code{NA} for other cells. Otherwise, the output values are \code{TRUE} and \code{FALSE}}
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
r <- rast(ncols=5, nrows=5, vals=1, ext=c(0,1,0,1))
r[10:20] <- NA
x <- not.na(r)
y <- not.na(r, falseNA=TRUE)
unique(values(c(x, y)))
}

\keyword{methods}
\keyword{spatial}

