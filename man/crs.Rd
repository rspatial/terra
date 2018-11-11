\name{crs}

\docType{methods}

\alias{crs}
\alias{crs<-}

\alias{crs,SpatRaster-method}
\alias{crs<-,SpatRaster-method}
\alias{crs<-,SpatRaster,character-method}

\alias{crs,SpatLayer-method}
\alias{crs<-,SpatLayer-method}
\alias{crs<-,SpatLayer,character-method}

\title{ Get or set a coordinate reference system (projection) }

\description{
Get or set the coordinate reference system (CRS) of a SpatRaster object. 
}

\usage{
\S4method{crs}{SpatRaster}(x)
\S4method{crs}{SpatLayer}(x)
\S4method{crs}{SpatRaster}(x, ...)<-value
\S4method{crs}{SpatLayer}(x, ...)<-value
}

\arguments{
  \item{x}{SpatRaster}
  \item{value}{character string describing a coordinate reference system in the PROJ.4 "+" format }
  \item{...}{additional arguments (none implemented)}
}


\value{
character or modified SpatRaster
}


\examples{ 
r <- rast()
crs(r)
#crs(r) <- "+proj=lcc +lat_1=48 +lat_2=33 +lon_0=-100 +ellps=WGS84"
#crs(r)
}

\keyword{ spatial }

