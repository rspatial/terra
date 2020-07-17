\name{crs}

\docType{methods}

\alias{crs}
\alias{crs<-}

\alias{crs,SpatRaster-method}
\alias{crs<-,SpatRaster-method}
\alias{crs<-,SpatRaster,ANY-method}

\alias{crs,SpatVector-method}
\alias{crs<-,SpatVector-method}
\alias{crs<-,SpatVector,ANY-method}

\title{ Get or set a coordinate reference system }

\description{
Get or set the coordinate reference system (CRS), also referred to as a "projection" of a SpatRaster or SpatVector object. 

Setting a new CRS does not change the data itself, it just changes the label. So you should only set the CRS of a dataset (if it does not come with one) to what it *is*, not to what you would *like it to be*. See \code{\link{project}} to *transform* spatial from one CRS to another.
}

\usage{
\S4method{crs}{SpatRaster}(x)

\S4method{crs}{SpatVector}(x)

\S4method{crs}{SpatRaster}(x, ...)<-value

\S4method{crs}{SpatVector}(x, ...)<-value
}

\arguments{
  \item{x}{SpatRaster or SpatVector}
  \item{value}{character string describing a coordinate reference system. This can be in a WKT format, as a EPSG code, or a PROJ.4 "+" format (but see Note)}
  \item{...}{additional arguments (none implemented)}
}


\note{
Because of changes in the PROJ library that is behind the coordinate transformations, when using the PROJ.4 format, the datum should be WGS84 -- if you want to transform your data to a different coordinate reference system)
}

\value{
character or modified SpatRaster/Vector
}


\examples{ 
r <- rast()
crs(r)
crs(r) <- "+proj=lcc +lat_1=48 +lat_2=33 +lon_0=-100 +ellps=WGS84"
crs(r)

# You can also use epsg codes
crs(r)  <- "epsg:25831"
}

\keyword{ spatial }

