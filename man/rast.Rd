\name{rast}

\docType{methods}

\alias{rast}
\alias{rast,missing-method}
\alias{rast,character-method}
\alias{rast,SpatRaster-method}
\alias{rast,matrix-method}
\alias{rast,array-method}
\alias{rast,SpatExtent-method}
\alias{rast,SpatVector-method}
\alias{rast,Raster-method}
\alias{rast,PackedSpatRaster-method}


\title{Create a SpatRaster}

\description{
Methods to create a SpatRaster. These objects can be created from scratch or from a file.

A SpatRaster represents a spatially referenced surface divided into three dimensional cells (rows, columns, and layers).

When a SpatRaster is created from a file, it does (initially) not contain any cell (pixel) values in memory (RAM), it only has the parameters that describe the SpatRaster. You can access cell-values with \code{\link{values}}.
}


\usage{
\S4method{rast}{character}(x, ...)

\S4method{rast}{missing}(x, nrows=180, ncols=360, nlyrs=1, xmin=-180, xmax=180, 
                                      ymin=-90, ymax=90, crs, extent, resolution, ...)

\S4method{rast}{SpatRaster}(x, nlyrs=nlyr(x), ...) 

\S4method{rast}{matrix}(x, crs=NA, type="", ...)

\S4method{rast}{SpatVector}(x, nrows=10, ncols=10, nlyrs=1, ...)
									  
\S4method{rast}{SpatExtent}(x, nrows=10, ncols=10, nlyrs=1, crs="", ...)

\S4method{rast}{Raster}(x, ...)
}

\arguments{
\item{x}{filename (character), Extent, SpatRaster, matrix, array, or Raster* object (from the "raster" package). Can also be missing}
\item{nrows}{positive integer. Number of rows}
\item{ncols}{positive integer. Number of columns}
\item{nlyrs}{positive integer. Number of layers}
\item{xmin}{minimum x coordinate (left border)}
\item{xmax}{maximum x coordinate (right border)}
\item{ymin}{minimum y coordinate (bottom border)}
\item{ymax}{maximum y coordinate (top border)}
\item{extent}{object of class SpatExtent. If present, the arguments xmin, xmax, ymin and ymax are ignored}
\item{crs}{character. PROJ.4 type description of a Coordinate Reference System (map projection). If this argument is missing, and the x coordinates are within -360 .. 360 and the y coordinates are within -90 .. 90, "+proj=longlat +datum=WGS84" is used}
\item{resolution}{numeric vector of length 1 or 2 to set the resolution (see \code{\link{res}}). If this argument is used, arguments \code{ncol} and \code{nrow} are ignored }
\item{type}{character. 
If the value is not \code{"xyz"}, the raster has the same number of rows and colums as the matrix.
If the value is \code{"xyz"}, the matrix must have at least two columns, the first with \code{x} (or longitude) and the second with \code{y} (or latitude) coordinates that represent the centers of raster cells. The additional columns are the values associated with the raster cells.}

\item{...}{additional arguments. None implemented }  
}


\value{
SpatRaster
}


\examples{
# Create a SpatRaster from scratch
x <- rast(nrow=108, ncol=21, xmin=0, xmax=10)

# Create a SpatRaster from a file
f <- system.file("exdata/test.tif", package="terra")
r <- rast(f)

s <- rast(system.file("exdata/logo.grd", package="terra"))   

# Create a skeleton with no associated cell values
rast(s)
}


\keyword{methods}
\keyword{spatial}

