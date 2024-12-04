\name{readAll}

\docType{methods}

\alias{readAll}
\alias{readAll,SpatRaster-method}
\alias{readAll,SpatRasterDataset-method}

\title{Read all values into memory}

\description{
Reads all values of a SpatRaster or SpatRasterDataset into memory.
}

\usage{
\S4method{readAll}{SpatRaster}(x)

\S4method{readAll}{SpatRasterDataset}(x)
}

\arguments{
  \item{x}{SpatRaster or SpatRasterDataset}
}

\value{
  Invisibily returns a logical value indicating whether the values were successfully read into memory.
}

\seealso{ \code{\link{set.values}} }


\keyword{ spatial }
\keyword{ methods }
