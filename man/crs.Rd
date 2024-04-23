\name{crs}

\docType{methods}

\alias{crs}
\alias{crs<-}

\alias{crs,SpatExtent-method}

\alias{crs,SpatRaster-method}
\alias{crs,SpatRasterDataset-method}
\alias{crs<-,SpatRaster-method}
\alias{crs<-,SpatRaster,ANY-method}

\alias{crs,SpatVector-method}
\alias{crs,SpatVectorProxy-method}
\alias{crs,SpatVectorCollection-method}
\alias{crs<-,SpatVector-method}
\alias{crs<-,SpatVector,ANY-method}

\alias{crs,sf-method}
\alias{crs,character-method}

\title{ Get or set a coordinate reference system }

\description{
Get or set the coordinate reference system (CRS), also referred to as a "projection", of a SpatRaster or SpatVector. 

Setting a new CRS does not change the data itself, it just changes the label. So you should only set the CRS of a dataset (if it does not come with one) to what it *is*, not to what you would *like it to be*. See \code{\link{project}} to *transform* an object from one CRS to another.
}

\usage{
\S4method{crs}{SpatRaster}(x, proj=FALSE, describe=FALSE, parse=FALSE)

\S4method{crs}{SpatVector}(x, proj=FALSE, describe=FALSE, parse=FALSE)

\S4method{crs}{character}(x, proj=FALSE, describe=FALSE, parse=FALSE)

\S4method{crs}{SpatRaster}(x, warn=FALSE)<-value

\S4method{crs}{SpatVector}(x, warn=FALSE)<-value
}

\arguments{
  \item{x}{SpatRaster or SpatVector}
  \item{proj}{logical. If \code{TRUE} the crs is returned in PROJ-string notation}
  \item{describe}{logical. If \code{TRUE} the name, EPSG code, and the name and extent of the area of use are returned if known}
  \item{warn}{logical. If \code{TRUE}, a message is printed when the object already has a non-empty crs}
  \item{value}{character string describing a coordinate reference system. This can be in a WKT format, as a <authority:number> code such as "EPSG:4326", or a PROJ-string format such as "+proj=utm +zone=12" (see Note)}
  \item{parse}{logical. If \code{TRUE}, wkt parts are parsed into a vector (each line becomes an element)}
}


\note{
Projections are handled by the PROJ/GDAL libraries. The PROJ developers suggest to define a CRS with the WKT2 or <authority>:<code> notation. It is not practical to define one's own custom CRS with WKT2, and the the <authority>:<code> system only covers a handful of (commonly used) CRSs. To work around this problem it is still possible to use the deprecated PROJ-string notation (\code{+proj=...}) with one major caveat: the datum should be WGS84 (or the equivalent NAD83) -- if you want to transform your data to a coordinate reference system with a different datum. Thus as long as you use WGS84, or an ellipsoid instead of a datum, you can safely use PROJ-strings to represent your CRS; including to define your own custom CRS.

You can also set the crs to "local" to get an informal coordinate system on an arbitrary Euclidean (Cartesian) plane with units in meter.
}


\value{
character or modified SpatRaster/Vector
}


\examples{ 
r <- rast()
crs(r)
crs(r, describe=TRUE, proj=TRUE)

crs(r) <- "+proj=lcc +lat_1=48 +lat_2=33 +lon_0=-100 +ellps=WGS84"
crs(r)

# You can use epsg codes
crs(r)  <- "epsg:25831"
crs(r, describe=TRUE)$area

crs("epsg:25831", describe=TRUE)
}

\keyword{ spatial }

