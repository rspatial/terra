\name{c}

\docType{methods}

\alias{c}
\alias{c,SpatRaster-method}
\alias{c,SpatStack-method}

\title{Combine}

\description{
Combine SpatRaster objects (similar to \code{raster::} \code{\link[raster]{stack}})
}

\usage{
\S4method{c}{SpatRaster}(x, ...)

\S4method{c}{SpatStack}(x, ...)
}

\arguments{
  \item{x}{SpatRaster or SpatStack}
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
