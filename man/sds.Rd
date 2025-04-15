\name{sds}

\docType{methods}

\alias{sds}
\alias{sds,missing-method}
\alias{sds,character-method}
\alias{sds,list-method}
\alias{sds,array-method}
\alias{sds,SpatRaster-method}
\alias{sds,stars-method}
\alias{sds,stars_proxy-method}
\alias{[<-,SpatRasterDataset,numeric,missing-method}


\title{Create a SpatRasterDataset}

\description{
Methods to create a SpatRasterDataset. This is an object to hold "sub-datasets", each represented by a SpatRaster that may have multiple layers. All sub-datasets must have the same raster geometry (extent and resolution). You can use a SpatRasterCollection (see \code{\link{sprc}}) to combine SpatRasters with different geometries.

See \code{\link{describe}} for getting information about the sub-datasets present in a file.
}

\usage{
\S4method{sds}{missing}(x) 

\S4method{sds}{character}(x, ids=0, opts=NULL, raw=FALSE, noflip=FALSE, guessCRS=TRUE, domains="")

\S4method{sds}{SpatRaster}(x, ...) 

\S4method{sds}{list}(x) 

\S4method{sds}{array}(x, crs="", extent=NULL)
}

\arguments{
\item{x}{character (filename), or SpatRaster, or list of SpatRasters, or missing. If multiple filenames are provided, it is attempted to make SpatRasters from these, and combine them into a SpatRasterDataset}
\item{ids}{optional. vector of integer subdataset ids. Ignored if the first value is not a positive integer}
\item{opts}{character. GDAL dataset open options}
\item{raw}{logical. If \code{TRUE}, scale and offset values are ignored}
\item{noflip}{logical. If \code{TRUE}, a raster (e.g. JPEG image) that is not georeferenced and that GDAL assigns a flipped extent to (\code{ymax < ymin}), is not considered flipped. This avoids the need to \code{\link{flip}} the raster vertically}
\item{guessCRS}{logical. If \code{TRUE} and the file does not specify a CRS but has an extent that is within longitude/latitude bounds, the longitude/latitude crs is assigned to the SpatRaster}
\item{domains}{character. Metadata domains to read (see \code{\link{metags}} to retrieve their values if there are any). "" is the default domain}
\item{crs}{character. Description of the Coordinate Reference System (map projection) in \code{PROJ.4}, \code{WKT} or \code{authority:code} notation. If this argument is missing, and the x coordinates are within -360 .. 360 and the y coordinates are within -90 .. 90, longitude/latitude is assigned}
\item{extent}{\code{\link{SpatExtent}}}
\item{...}{additional \code{SpatRaster} objects} 
}


\value{
SpatRasterDataset
}

\seealso{
\code{\link{sprc}}, \code{\link{describe}}
}

\examples{
s <- rast(system.file("ex/logo.tif", package="terra"))   
x <- sds(s, s/2)
names(x) <- c("first", "second")
x
length(x)

# extract the second SpatRaster
x[2]

a <- array(1:9, c(3,3,3,3))
sds(a)
}


\keyword{methods}
\keyword{spatial}

