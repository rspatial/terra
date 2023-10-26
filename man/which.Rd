\name{which.lyr}

\docType{methods}

\alias{which.lyr}
\alias{which.lyr,SpatRaster-method}

\title{Which cells are TRUE?}

\description{
This method returns a single layer SpatRaster with cell values indicating the first layer in the input that is \code{TRUE}. All numbers that are not zero (or \code{FALSE}), are considered to be \code{TRUE}.
}


\usage{
\S4method{which.lyr}{SpatRaster}(x)
}

\arguments{
  \item{x}{SpatRaster}
}


\value{
SpatRaster
}

\seealso{\code{\link{isTRUE}}, \code{\link{which}}, See \code{\link{Summary-methods}} for \code{which.min} and \code{which.max}}


\examples{
s <- rast(system.file("ex/logo.tif", package="terra"))   
x <- which.lyr(s > 100)
}

\keyword{spatial}
