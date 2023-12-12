\name{rast}

\docType{methods}

\alias{rast}
\alias{rast,missing-method}
\alias{rast,character-method}
\alias{rast,SpatRaster-method}
\alias{rast,SpatRasterDataset-method}
\alias{rast,data.frame-method}
\alias{rast,matrix-method}
\alias{rast,array-method}
\alias{rast,list-method}
\alias{rast,SpatExtent-method}
\alias{rast,SpatVector-method}
\alias{rast,ANY-method}
\alias{rast,PackedSpatRaster-method}
\alias{rast,stars-method}
\alias{rast,stars_proxy-method}

\title{Create a SpatRaster}

\description{
Methods to create a SpatRaster. These objects can be created from scratch, from a filename, or from another object.

A SpatRaster represents a spatially referenced surface divided into three dimensional cells (rows, columns, and layers). 

When a SpatRaster is created from a file, it does not load the cell (pixel) values into memory (RAM). It only reads the parameters that describe the geometry of the SpatRaster, such as the number of rows and columns and the coordinate reference system. The actual values will be read when needed.
}

\usage{
\S4method{rast}{character}(x, subds=0, lyrs=NULL, drivers=NULL, opts=NULL, 
		  win=NULL, snap="near", vsi=FALSE, raw=FALSE)

\S4method{rast}{missing}(x, nrows=180, ncols=360, nlyrs=1, xmin=-180, xmax=180, 
          ymin=-90, ymax=90, crs, extent, resolution, vals, names, time, units)

\S4method{rast}{SpatRaster}(x, nlyrs=nlyr(x), names, vals, keeptime=TRUE, 
		  keepunits=FALSE, props=FALSE, tags=FALSE) 

\S4method{rast}{matrix}(x, type="", crs="", digits=6, extent=NULL)

\S4method{rast}{data.frame}(x, type="xyz", crs="", digits=6, extent=NULL)

\S4method{rast}{array}(x, crs="", extent=NULL)

\S4method{rast}{list}(x, warn=TRUE)

\S4method{rast}{SpatRasterDataset}(x)

\S4method{rast}{SpatVector}(x, ...)
									  
\S4method{rast}{SpatExtent}(x, ...)
}

\arguments{
\item{x}{filename (character), missing, SpatRaster, SpatRasterDataset, SpatExtent, SpatVector, matrix, array, list of SpatRaster objects. For other types it will be attempted to create a SpatRaster via (`as(x, "SpatRaster")`}
\item{subds}{positive integer or character to select a sub-dataset. If zero or "", all sub-datasets are returned (if possible)}
\item{lyrs}{positive integer or character to select a subset of layers (a.k.a. "bands")}
\item{drivers}{character. GDAL drivers to consider}
\item{opts}{character. GDAL dataset open options}
\item{win}{SpatExtent to set a \code{\link{window}} (area of interest)}
\item{snap}{character. One of "near", "in", or "out", to indicate how the extent of \code{\link{window}} should be "snapped" to \code{x}}
\item{vsi}{logical. If \code{TRUE}, "\\vsicurl\\" is prepended to filenames that start with "http"}
\item{raw}{logical. If \code{TRUE}, scale and offset values are ignored}
\item{nrows}{positive integer. Number of rows}
\item{ncols}{positive integer. Number of columns}
\item{nlyrs}{positive integer. Number of layers}
\item{xmin}{minimum x coordinate (left border)}
\item{xmax}{maximum x coordinate (right border)}
\item{ymin}{minimum y coordinate (bottom border)}
\item{ymax}{maximum y coordinate (top border)}
\item{crs}{character. Description of the Coordinate Reference System (map projection) in \code{PROJ.4}, \code{WKT} or \code{authority:code} notation. See \code{\link{crs}}. If this argument is missing, and the x coordinates are within -360 .. 360 and the y coordinates are within -90 .. 90, longitude/latitude is assigned}
\item{keeptime}{logical. If \code{FALSE} the time stamps are discarded}
\item{keepunits}{logical. If \code{FALSE} the layer units are discarded}
\item{props}{logical. If \code{TRUE} the properties (categories and color-table) are kept}
\item{tags}{logical. If \code{TRUE} the user specified metadata tags are kept (see \code{\link{metags}}).}
\item{extent}{object of class SpatExtent. If present, the arguments xmin, xmax, ymin and ymax are ignored}
\item{resolution}{numeric vector of length 1 or 2 to set the spatial resolution (see \code{\link{res}}). If this argument is used, arguments \code{ncols} and \code{nrows} are ignored }
\item{vals}{numeric. An optional vector with cell values (if fewer values are provided, these are recycled to reach the number of cells)}
\item{names}{character. An optional vector with layer names (must match the number of layers)}
\item{time}{time or date stamps for each layer}
\item{units}{character. units for each layer}
\item{type}{character. If the value is \code{"xyz"}, the matrix or data.frame \code{x} must have at least two columns, the first with \code{x} (or longitude) and the second with \code{y} (or latitude) coordinates that represent the centers of raster cells. The additional columns are the values associated with the raster cells. If the value is \code{"xylz"}, \code{x} must have four columns with the third representing the layer and the fourth the cell values. If the value is \code{""}, the resulting SpatRaster will have the same number of rows and columns as \code{x}. 
}
\item{digits}{integer to set the precision for detecting whether points are on a regular grid (a low number of digits is a low precision). Only used when \code{type="xyz"}}
\item{warn}{logical. If \code{TRUE}, a warnings about empty rasters may be emitted}
\item{...}{additional arguments passed on to the \code{rast,missing-method}}  
}


\value{
SpatRaster
}

\details{
Files are read with the GDAL library. GDAL guesses the file format from the name, and/or tries reading it with different "drivers" (see \code{\link{gdal}}) until it succeeds. In very few cases this may cause a file to be opened with the wrong driver, and some information may be lost. For example, when a netCDF file is opened with the HDF5 driver. You can avoid that by using argument \code{rast("filename.ncdf", drivers="NETCDF")}

These classes hold a C++ pointer to the data "reference class" and that creates some limitations. They cannot be recovered from a saved R session either or directly passed to nodes on a computer cluster. Generally, you should use \code{\link{writeRaster}} to save SpatRaster objects to disk (and pass a filename or cell values ot cluster nodes). Also see \code{\link{wrap}}. 
}

\seealso{\code{\link{sds}} to create a SpatRasterDataset (4 dimensions) and \code{\link{vect}} for vector (points, lines, polygons) data}

\examples{
# Create a SpatRaster from scratch
x <- rast(nrows=108, ncols=21, xmin=0, xmax=10)

# Create a SpatRaster from a file
f <- system.file("ex/elev.tif", package="terra")
r <- rast(f)

s <- rast(system.file("ex/logo.tif", package="terra"))   

# Create a skeleton with no associated cell values
rast(s)

# from a matrix 
m <- matrix(1:25, nrow=5, ncol=5)
rm <- rast(m)

# from a "xyz" data.frame
d <- as.data.frame(rm, xy=TRUE)
head(d)
rast(d, type="xyz")

}

\keyword{methods}
\keyword{spatial}

