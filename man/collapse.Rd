\name{collapse}

\docType{methods}

\alias{collapse}
\alias{collapse,SpatRaster-method}
\alias{collapse,SpatRasterDataset-method}

\title{Collapse SpatRaster or SpatRasterDataset objects}

\description{
Combines sources within a SpatRaster object (that are in memory, or from the same file) to allow for faster processing. 

Or combine subdatsets into a SpatRaster.
}
\usage{

\S4method{collapse}{SpatRaster}(x, ...)

\S4method{collapse}{SpatRasterDataset}(x)
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
collapse(x)
}

\keyword{spatial}
