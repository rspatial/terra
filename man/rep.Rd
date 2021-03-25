\name{rep}

\docType{methods}

\alias{rep}
\alias{rep,SpatRaster-method}

\title{Replicate layers}

\description{
Replicate layers in a SpatRaster
}

\usage{
\S4method{rep}{SpatRaster}(x, ...)

}

\arguments{
  \item{x}{SpatRaster}
  \item{...}{arguments as in \code{\link[base]{rep}}}
}



\value{
SpatRaster
}


\examples{
s <- rast(system.file("ex/logo.tif", package="terra"))   
x <- rep(s, 2)
nlyr(x)
names(x)
x
}

\keyword{spatial}
