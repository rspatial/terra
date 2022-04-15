\name{where}

\docType{methods}

\alias{where.min}
\alias{where.min,SpatRaster-method}
\alias{where.max}
\alias{where.max,SpatRaster-method}

\title{Where are the cells with the min or max values?}

\description{
This method returns the cell numbers for the cells with the min or max values of a SpatRaster.
}


\usage{
\S4method{where.min}{SpatRaster}(x)
\S4method{where.max}{SpatRaster}(x)
}

\arguments{
  \item{x}{SpatRaster}
}


\value{
list
}

\seealso{\code{\link{which}} and \code{\link{Summary-methods}} for \code{which.min} and \code{which.max}}


\examples{
f <- system.file("ex/elev.tif", package="terra")
r <- rast(f)
x <- where.min(r)
}

\keyword{spatial}
