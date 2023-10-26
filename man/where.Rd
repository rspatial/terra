\name{where}

\docType{methods}

\alias{where.min}
\alias{where.min,SpatRaster-method}
\alias{where.max}
\alias{where.max,SpatRaster-method}

\title{Where are the cells with the min or max values?}

\description{
This method returns the cell numbers for the cells with the min or max values of each layer in a SpatRaster.
}


\usage{
\S4method{where.min}{SpatRaster}(x, values=TRUE, list=FALSE)

\S4method{where.max}{SpatRaster}(x, values=TRUE, list=FALSE)
}

\arguments{
  \item{x}{SpatRaster}
  \item{values}{logical. If \code{TRUE} the min or max values are also returned}
  \item{list}{logical. If \code{TRUE} a list is returned instead of a matrix}
}


\value{
matrix or list
}

\seealso{\code{\link{which}} and \code{\link{Summary-methods}} for \code{which.min} and \code{which.max}}


\examples{
f <- system.file("ex/elev.tif", package="terra")
r <- rast(f)
where.min(r)
}

\keyword{spatial}
