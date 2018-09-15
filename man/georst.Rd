\name{georst}

\docType{methods}

\alias{georst}
\alias{georst,missing-method}
\alias{georst,character-method}
\alias{georst,SpatRaster-method}
\alias{georst,matrix-method}
\alias{georst,array-method}

\title{Create a SpatRaster object}

\description{
Methods to create a SpatRaster object. These objects can be created from scratch or from a file.

A SpatRaster represents a spatially referenced surface divided into three dimensional cells (rows, columns, and layers).

When a SpatRaster is created from a file, it does (initially) not contain any cell (pixel) values in (RAM) memory, it only has the parameters that describe the SpatRaster. You can access cell-values with \code{\link{values}}.
}


\usage{

\S4method{georst}{missing}(nrow=180, ncol=360, nlayer=1, xmin=-180, xmax=180, ymin=-90, ymax=90, crs, extent, resolution, ...)

\S4method{georst}{character}(x, ...)

\S4method{georst}{SpatRaster}(x, ...) 
}

\arguments{
\item{x}{filename (character), Extent, SpatRaster, matrix, or array}
\item{...}{Additional arguments, see Details }  
\item{nrow}{positive integer. Number of rows}
\item{ncol}{postive integer. Number of columns}
\item{nlayer}{postive integer. Number of layers}
\item{xmin}{minimum x coordinate (left border)}
\item{xmax}{maximum x coordinate (right border)}
\item{ymin}{minimum y coordinate (bottom border)}
\item{ymax}{maximum y coordinate (top border)}
\item{extent}{object of class SpatExtent. If present, the arguments xmn, xmx, ymn and ynx are ignored}
\item{crs}{character. PROJ.4 type description of a Coordinate Reference System (map projection). If this argument is missing, and the x coordinates are within -360 .. 360 and the y coordinates are within -90 .. 90, "+proj=longlat +datum=WGS84" is used. Also see under Details if \code{x} is a character (filename)}
\item{resolution}{numeric vector of length 1 or 2 to set the resolution (see \code{\link{res}}). If this argument is used, arguments \code{ncol} and \code{nrow} are ignored }
}


\value{
SpatRaster
}


\examples{
# Create a SpatRaster object from scratch
r1 <- georst(nrow=108, ncol=21, xmin=0, xmax=10)

# Create a SpatRaster object from a file
#   N.B.: For your own files, omit the 'system.file' and 'package="terra"' bits
#   these are just to get the path to files installed with the package

f <- system.file("external/test.grd", package="terra")
f
r <- georst(f)

logo <- georst(system.file("external/rlogo.grd", package="raster"))   
}


\keyword{methods}
\keyword{spatial}

