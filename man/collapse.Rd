\name{tighten}

\docType{methods}

\alias{tighten}
\alias{tighten,SpatRaster-method}
\alias{tighten,SpatRasterDataset-method}

\title{tighten SpatRaster or SpatRasterDataset objects}

\description{
Combines data sources within a SpatRaster (that are in memory, or from the same file) to allow for faster processing. 

Or combine sub-datasets into a SpatRaster.
}
\usage{

\S4method{tighten}{SpatRaster}(x)

\S4method{tighten}{SpatRasterDataset}(x)
}

\arguments{
  \item{x}{SpatRaster or SpatRasterDataset}
}


\value{
SpatRaster
}

\examples{
r <- rast(nrow=5, ncol=9, vals=1:45)
x <- c(r, r*2, r*3)
x
tighten(x)
}

\keyword{spatial}
