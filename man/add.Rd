\name{add}

\docType{methods}

\alias{add<-}
\alias{add<-,SpatRaster,SpatRaster-method}
\alias{add<-,SpatRasterCollection,SpatRaster-method}
\alias{add<-,SpatRasterDataset,SpatRaster-method}


\title{Add (in place) a SpatRaster to another SpatRaster object or to a SpatRasterDataset or SpatRasterCollection}

\description{
Add (in place) a SpatRaster to another SpatRaster object. Comparable with \code{\link{c}}, but without copying the object. 
}

\usage{
\S4method{add}{SpatRaster,SpatRaster}(x)<-value

\S4method{add}{SpatRasterDataset,SpatRaster}(x)<-value

\S4method{add}{SpatRasterCollection,SpatRaster}(x)<-value
}

\arguments{
  \item{x}{SpatRaster, SpatRasterDataset or SpatRasterCollection}
  \item{value}{SpatRaster}
}

\seealso{\code{\link{c}}}

\value{
SpatRaster
}

\examples{
r <- rast(nrows=5, ncols=9, vals=1:45)
x <- c(r, r*2)
add(x) <- r*3
x
}

\keyword{spatial}
