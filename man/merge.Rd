\name{merge}

\docType{methods}

\alias{merge}
\alias{merge,SpatRaster,SpatRaster-method}
\alias{merge,SpatRasterCollection,missing-method}
\alias{merge,SpatVector,SpatVector-method}
\alias{merge,SpatVector,data.frame-method}

\title{ 
Merge SpatRasters, or merge a SpatVector with a data.frame
}

\description{
Merge multiple SpatRasters to create a new SpatRaster with a larger spatial extent. The SpatRasters should all have the same coordinate reference system. They should normally also have the same spatial origin and resolution, but automatic resampling can be done depending on the algorithm used (see argument \code{algo}). In areas where the SpatRasters overlap, the values of the SpatRaster that is first in the sequence of arguments (or in the SpatRasterCollection) will be retained (unless \code{first=FALSE}). 

There is also a method for merging SpatVector with a data.frame; that is, to join the data.frame to the attribute table of the SpatVector.
}

\usage{
\S4method{merge}{SpatRaster,SpatRaster}(x, y, ..., first=TRUE, na.rm=TRUE, algo=1, method=NULL, 
			filename="", overwrite=FALSE, wopt=list())

\S4method{merge}{SpatRasterCollection,missing}(x, first=TRUE, na.rm=TRUE, algo=1, method=NULL, filename="", ...)

\S4method{merge}{SpatVector,data.frame}(x, y, ...)
}

\arguments{
  \item{x}{SpatRaster, SpatRasterCollection, or SpatVector}
  \item{y}{missing if \code{x} is a SpatRasterCollection. SpatRaster if \code{x} is a SpatRaster. data.frame if \code{x} is a SpatVector}
  \item{...}{if \code{x} is a SpatRaster: additional objects of the same class as \code{x}. If \code{x} is a SpatRasterCollection: options for writing files as in \code{\link{writeRaster}}. If \code{x} is a SpatVector, the same arguments as in \code{\link[base]{merge}}}
  \item{first}{logical. If \code{TRUE}, in areas where rasters overlap, the first value is used. Otherwise the last value is used}
  \item{na.rm}{logical. If \code{TRUE} missing values are are ignored. This is only used for algo 1; the other two always ignore missing values}
  \item{algo}{integer. You can use 1, 2 or 3 to pick a merge algorithm. algo 1 is generally faster than algo 2, but it may have poorer file compression. Algo 1 will resample input rasters (and that may slow it down), but algo 2 does not do that. You can increase the tolerance option to effectively get nearest neighbor resampling with, for example, \code{wopt=list(tolerance=.2)} allows misalignment of .2 times the resolution of the first input raster and effectively use nearest neighbor resampling. Algo 3 creates a virtual raster (see \code{\link{vrt}}) and returned if no filename is specified. This is very quick and can be a good approach if the merge raster is used as input to a next step in the analysis. It allows any amount of misalignment (and does not respond to the tolerance option). Otherwise its speed is like that of algo 2}
  \item{method}{character. The interpolation method to be used if resampling is necessary (see argument \code{algo}). One of "nearest", "bilinear", "cubic", "cubicspline", "lanczos", "average", "mode" as in \code{\link{resample}}. If \code{NULL}, "nearest" is used for categorical rasters and "bilinear" for other rasters} 
  \item{filename}{character. Output filename}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list with named options for writing files as in \code{\link{writeRaster}}}
} 


\value{
SpatRaster or SpatVector
}

\seealso{
Combining tiles with \code{\link{vrt}} may be more efficient than using \code{merge}. See \code{\link{mosaic}} for averaging overlapping regions. 

See \code{\link{classify}} to merge a \code{SpatRaster} and a \code{data.frame} and \code{\link{union}} to combine SpatExtent objects. 
}


\examples{
x <- rast(xmin=-110, xmax=-80, ymin=40, ymax=70, res=1, vals=1)
y <- rast(xmin=-85, xmax=-55, ymax=60, ymin=30, res=1, vals=2)
z <- rast(xmin=-60, xmax=-30, ymax=50, ymin=20, res=1, vals=3)

m1 <- merge(x, y, z)
m2 <- merge(z, y, x)
m3 <- merge(y, x, z)
# panel(c(m1, m2, m3))

# if you have many SpatRasters, it may be convenient
# to make a SpatRasterCollection
# s <- sprc(list(x, y, z))
s <- sprc(x, y, z)

sm1 <- merge(s, algo=1, first=FALSE)
sm2 <- merge(s, algo=2, first=FALSE)
#sm3 <- merge(s, algo=3, first=FALSE)

## SpatVector with data.frame
f <- system.file("ex/lux.shp", package="terra")
p <- vect(f)
dfr <- data.frame(District=p$NAME_1, Canton=p$NAME_2, Value=round(runif(length(p), 100, 1000)))
dfr <- dfr[1:5, ]
pm <- merge(p, dfr, all.x=TRUE, by.x=c('NAME_1', 'NAME_2'), by.y=c('District', 'Canton'))
pm
values(pm)
}

\keyword{methods}
\keyword{spatial}

