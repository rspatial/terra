\name{extent}

\docType{methods}


\alias{ext}
\alias{ext<-}
\alias{ext,GeoExtent-method}
\alias{ext,SpatRaster-method}
\alias{ext<-,SpatRaster,GeoExtent-method}
\alias{ext,numeric-method}
\alias{ext,missing-method}

  
\title{ext}

\description{
Get a GeoExtent object of a SpatRaster object, or coordinates from such an object. Or create a GeoExtent object from a vector (length=4; order= xmin, xmax, ymin, ymax)
}

\usage{
\S4method{ext}{SpatRaster}(x, ...)
\S4method{ext}{SpatRaster,GeoExtent}(x)<-value
}

\arguments{
  \item{x}{SpatRaster object}
  \item{value}{GeoExtent object}  
  \item{...}{Additional arguments. None implemented }
}
 
\value{
GeoExtent object or numeric coordinate
}
  
\author{Robert J. Hijmans}


\examples{
r <- georst()
ext(r)
ext(c(0, 20, 0, 20))

}

\keyword{spatial}

