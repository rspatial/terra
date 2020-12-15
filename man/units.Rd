\name{units}

\alias{units}
\alias{units<-}
\alias{units,SpatRaster-method}
\alias{units<-,SpatRaster-method}
\alias{units,SpatRasterDataset-method}
\alias{units<-,SpatRasterDataset-method}

  
\title{units of SpatRaster or SpatRasterDataSet}

\description{
Get or set the units of the layers of a SpatRaster or the datasets in a SpatRasterDataSet.
}

\usage{
\S4method{units}{SpatRaster}(x)

\S4method{units}{SpatRaster}(x)<-value

\S4method{units}{SpatRasterDataset}(x)

\S4method{units}{SpatRasterDataset}(x)<-value
}

\arguments{
  \item{x}{SpatRaster}
  \item{value}{character}
}

\value{
character  
}

\seealso{\code{\link{time}, \link{names}}}

\examples{
s <- rast(system.file("ex/logo.tif", package="terra"))   

units(s) <- c("m/s", "kg", "ha")
units(s)
s

units(s) <- "kg"
units(s)
}


\keyword{spatial}

