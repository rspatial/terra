\name{which.lyr}

\docType{methods}

\alias{which.lyr}
\alias{which.lyr,SpatRaster-method}

\title{Which cells are TRUE?}

\description{
This method returns a single layer SpatRaster with cell values that are the first layer in the input that has the value is not zero (FALSE), and, hence, is \code{TRUE}.
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
