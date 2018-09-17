\name{extent}

\docType{methods}


\alias{ext}
\alias{ext<-}
\alias{ext,SpatExtent-method}
\alias{ext,SpatRaster-method}
\alias{ext,SpatVector-method}
\alias{ext<-,SpatRaster,SpatExtent-method}
\alias{ext,numeric-method}
\alias{ext,missing-method}

  
\title{ext}

\description{
Get a SpatExtent object of a SpatRaster object, or coordinates from such an object. Or create a SpatExtent object from a vector (length=4; order= xmin, xmax, ymin, ymax)
}

\usage{
\S4method{ext}{SpatRaster}(x, ...)
\S4method{ext}{SpatRaster,SpatExtent}(x)<-value
}

\arguments{
  \item{x}{SpatRaster object}
  \item{value}{SpatExtent object}  
  \item{...}{Additional arguments. None implemented }
}
 
\value{
SpatExtent object or numeric coordinate
}
  
\author{Robert J. Hijmans}


\examples{
r <- rast()
ext(r)
ext(c(0, 20, 0, 20))

}

\keyword{spatial}

