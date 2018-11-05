\name{rast}

\docType{methods}

\alias{rast}
\alias{rast,missing-method}
\alias{rast,character-method}
\alias{rast,SpatRaster-method}
\alias{rast,matrix-method}
\alias{rast,array-method}
\alias{rast,SpatExtent-method}
\alias{rast,SpatLayer-method}


\title{Create a SpatRaster object}

\description{
Methods to create a SpatRaster object. These objects can be created from scratch or from a file.

A SpatRaster represents a spatially referenced surface divided into three dimensional cells (rows, columns, and layers).

When a SpatRaster is created from a file, it does (initially) not contain any cell (pixel) values in memory (RAM), it only has the parameters that describe the SpatRaster. You can access cell-values with \code{\link{values}}.
}


\usage{
\S4method{rast}{character}(x, ...)

\S4method{rast}{missing}(x, nrow=180, ncol=360, nlyr=1, xmin=-180, xmax=180, 
                                      ymin=-90, ymax=90, crs, extent, resolution, ...)

\S4method{rast}{SpatRaster}(x, ...) 

\S4method{rast}{SpatLayer}(x, nrow=10, ncol=10, nlyr=1, ...)
									  
\S4method{rast}{SpatExtent}(x, nrow=10, ncol=10, nlyr=1, crs="", ...)
}

\arguments{
\item{x}{filename (character), Extent, SpatRaster, matrix, or array. Can also be missing}
\item{nrow}{positive integer. Number of rows}
\item{ncol}{positive integer. Number of columns}
\item{nlyr}{positive integer. Number of layers}
\item{xmin}{minimum x coordinate (left border)}
\item{xmax}{maximum x coordinate (right border)}
\item{ymin}{minimum y coordinate (bottom border)}
\item{ymax}{maximum y coordinate (top border)}
\item{extent}{object of class SpatExtent. If present, the arguments xmin, xmax, ymin and ymax are ignored}
\item{crs}{character. PROJ.4 type description of a Coordinate Reference System (map projection). If this argument is missing, and the x coordinates are within -360 .. 360 and the y coordinates are within -90 .. 90, "+proj=longlat +datum=WGS84" is used}
\item{resolution}{numeric vector of length 1 or 2 to set the resolution (see \code{\link{res}}). If this argument is used, arguments \code{ncol} and \code{nrow} are ignored }
\item{...}{Additional arguments. None implemented }  
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

