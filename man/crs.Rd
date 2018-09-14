\name{crs}

\docType{methods}

\alias{crs}
\alias{crs<-}
\alias{crs,SpatRaster-method}
\alias{crs<-,SpatRaster-method}
\alias{crs<-,SpatRaster,character-method}


\title{ Get or set a coordinate reference system (projection) }

\description{
Get or set the coordinate reference system (CRS) of a SpatRaster object. 
}

\usage{
\S4method{crs}{SpatRaster}(x)

\S4method{crs}{SpatRaster}(x, ...)<-value
}

\arguments{
  \item{x}{SpatRaster object }
  \item{...}{additional arguments (not implemented)}
  \item{value}{character string describing a coordinate reference system in the PROJ.4 "+" format }
}


\value{
character or modified SpatRaster
}


\examples{ 
r <- georst()
crs(r)
#crs(r) <- "+proj=lcc +lat_1=48 +lat_2=33 +lon_0=-100 +ellps=WGS84"
#crs(r)
}

\keyword{ spatial }

