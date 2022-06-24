\name{sort}

\docType{methods}

\alias{sort}
\alias{sort,SpatRaster-method}

\title{Sort a SpatRaster}

\description{
Sort the cell values of a SpatRaster across layers
}

\usage{
\S4method{sort}{SpatRaster}(x, decreasing=FALSE, filename="", ...)

}

\arguments{
  \item{x}{SpatRaster}
  \item{decreasing}{logical}
  \item{filename}{character. Output filename}
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
}


\value{
SpatRaster
}


\examples{
f <- system.file("ex/elev.tif", package="terra")
r <- rast(f)
r <- c(r, r/2, r*2)
sort(r)
}

\keyword{spatial}

