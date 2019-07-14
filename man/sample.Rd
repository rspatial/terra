\name{sampleRegular}

\docType{methods}

\alias{sampleRegular}
\alias{sampleRegular,SpatRaster,numeric-method}

\title{Take a regular sample}

\description{
Take a regular sample of a SpatRaster. The returned object has the same extent, but fewer cells. 
}

\usage{
\S4method{sampleRegular}{SpatRaster,numeric}(x, size, ...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{size}{numeric. The sample size}
  \item{...}{additional arguments. None implemented}
}



\value{
SpatRaster
}


\examples{
f <- system.file("exdata/test.tif", package="terra")
r <- rast(f)
s <- sampleRegular(r, 1000)
}

\keyword{spatial}
