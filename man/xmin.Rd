\name{xmin}

\docType{methods}

\alias{xmin}
\alias{xmax}
\alias{ymin}
\alias{ymax}
\alias{xmin<-}
\alias{xmax<-}
\alias{ymin<-}
\alias{ymax<-}
\alias{xmin,SpatExtent-method}
\alias{xmin,SpatRaster-method}
\alias{xmin,SpatVector-method}
\alias{xmax,SpatExtent-method}
\alias{xmax,SpatRaster-method}
\alias{xmax,SpatVector-method}
\alias{ymin,SpatExtent-method}
\alias{ymin,SpatRaster-method}
\alias{ymin,SpatVector-method}
\alias{ymax,SpatExtent-method}
\alias{ymax,SpatRaster-method}
\alias{ymax,SpatVector-method}

\alias{xmin<-,SpatExtent,numeric-method}
\alias{xmin<-,SpatRaster,numeric-method}
\alias{xmax<-,SpatExtent,numeric-method}
\alias{xmax<-,SpatRaster,numeric-method}
\alias{ymin<-,SpatExtent,numeric-method}
\alias{ymin<-,SpatRaster,numeric-method}
\alias{ymax<-,SpatExtent,numeric-method}
\alias{ymax<-,SpatRaster,numeric-method}

  
\title{Get or set single values of an extent}

\description{
Get or set single values of an extent. Values can be set for a SpatExtent or SpatRaster, but not for a SpatVector)
}

\usage{
\S4method{xmin}{SpatExtent}(x)

\S4method{xmax}{SpatExtent}(x)

\S4method{ymin}{SpatExtent}(x)

\S4method{ymax}{SpatExtent}(x)

\S4method{xmin}{SpatRaster}(x)

\S4method{xmax}{SpatRaster}(x)

\S4method{ymin}{SpatRaster}(x)

\S4method{ymax}{SpatRaster}(x)

\S4method{xmin}{SpatVector}(x)

\S4method{xmax}{SpatVector}(x)

\S4method{ymin}{SpatVector}(x)

\S4method{ymax}{SpatVector}(x)

\S4method{xmin}{SpatRaster,numeric}(x)<-value

\S4method{xmax}{SpatRaster,numeric}(x)<-value

\S4method{ymin}{SpatRaster,numeric}(x)<-value

\S4method{ymax}{SpatRaster,numeric}(x)<-value
}

\arguments{
  \item{x}{SpatRaster, SpatExtent, or SpatVector}
  \item{value}{numeric}  
}
 
\value{
SpatExtent or numeric coordinate
}
  


\examples{
r <- rast()
ext(r)
ext(c(0, 20, 0, 20))

xmin(r)
xmin(r) <- 0
xmin(r)
}

\keyword{spatial}

