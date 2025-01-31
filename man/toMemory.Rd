\name{toMemory}

\docType{methods}

\alias{toMemory}
\alias{toMemory,SpatRaster-method}
\alias{toMemory,SpatRasterDataset-method}

\title{Read all cell values into memory}

\description{
Reads all cell values of a SpatRaster or SpatRasterDataset into memory.

Using this method is discouraged as it is not necessary for processing the data and may lead to excessive memory use that will slow down your computer or worse. It cannot be used for SpatRasters that are based on very large files.

The method may be useful if a relatively small dataset is used repeatedly, such that efficiency gains are made because the values only need to be read from disk once.
}

\usage{
\S4method{toMemory}{SpatRaster}(x)

\S4method{toMemory}{SpatRasterDataset}(x)
}

\arguments{
  \item{x}{SpatRaster or SpatRasterDataset}
}

\value{
Same as \code{x}
}


\seealso{ \code{\link{values}, \link{as.data.frame}, \link{readValues}, \link{inMemory}} }


\examples{
f <- system.file("ex/elev.tif", package="terra")
r <- rast(f)
sources(r)
inMemory(r)
x <- toMemory(r)
inMemory(x)
}

\keyword{ spatial }
\keyword{ methods }
