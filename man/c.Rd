\name{c}

\docType{methods}

\alias{c}
\alias{c, SpatRaster-method}

\title{Combine}

\description{
Combine SpatRaster objects (similar to raster::stack).
}

\usage{
\S4method{c}{SpatRaster}(x, ...)

}

\arguments{
  \item{x}{SpatRaster}
  \item{...}{SpatRaster objects to be combined with \code{x}}
}



\value{
SpatRaster
}


\examples{
r <- rast(nrow=5, ncol=9)
values(r) <- 1:ncell(r)
x <- c(r, r*2, r*3)
}

\keyword{spatial}
