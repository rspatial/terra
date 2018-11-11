\name{extent}

\docType{methods}


\alias{ext}
\alias{ext<-}
\alias{ext,SpatExtent-method}
\alias{ext,SpatRaster-method}
\alias{ext,SpatLayer-method}
\alias{ext<-,SpatRaster,SpatExtent-method}
\alias{ext,numeric-method}
\alias{ext,missing-method}
\alias{as.character,SpatExtent-method}

  
\title{ext}

\description{
Get a SpatExtent of a SpatRaster, or coordinates from such an object. Or create a SpatExtent from a vector (length=4; order= xmin, xmax, ymin, ymax)
}

\usage{
\S4method{ext}{SpatRaster}(x, ...)
\S4method{ext}{SpatRaster,SpatExtent}(x)<-value
\S4method{as.character}{SpatExtent}(x, ...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{value}{SpatExtent}  
  \item{...}{Additional arguments. None implemented }
}
 
\value{
SpatExtent or numeric coordinate
}
  
\author{Robert J. Hijmans}


\examples{
r <- rast()
ext(r)
ext(c(0, 20, 0, 20))

}

\keyword{spatial}

