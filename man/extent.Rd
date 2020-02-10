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
\alias{bbox,SpatRaster-method}
\alias{bbox,SpatVector-method}


  
\title{Create, get or set a SpatExtent}

\description{
Get a SpatExtent of a SpatRaster, or coordinates from such an object. Or create a SpatExtent from a vector (length=4; order= xmin, xmax, ymin, ymax)
}

\usage{
\S4method{ext}{SpatRaster}(x, ...)
\S4method{ext}{SpatRaster,SpatExtent}(x)<-value
}

\arguments{
  \item{x}{SpatRaster}
  \item{value}{SpatExtent}  
  \item{...}{additional arguments. None implemented }
}
 
\value{
SpatExtent or numeric coordinate
}
  


\examples{
r <- rast()
as.character(r)
as.character(ext(r))
}

\keyword{spatial}

