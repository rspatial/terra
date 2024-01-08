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
Methods to create a SpatRasterDataset. This is an object to hold "sub-datasets", each a SpatRaster that in most cases will have multiple layers. 

See \code{\link{describe}} for getting information about the sub-datasets present in a file.
}

\usage{
\S4method{sds}{missing}(x) 

\S4method{sds}{character}(x, ids=0, opts=NULL, raw=FALSE)

\S4method{sds}{SpatRaster}(x, ...) 

\S4method{sds}{list}(x) 

\S4method{sds}{array}(x, crs="", extent=NULL)
}

\arguments{
\item{x}{character (filename), or SpatRaster, or list of SpatRasters, or missing. If multiple filenames are provided, it is attempted to make SpatRasters from these, and combine them into a SpatRasterDataset}
\item{ids}{optional. vector of integer subdataset ids. Ignored if the first value is not a positive integer}
\item{opts}{character. GDAL dataset open options}
\item{raw}{logical. If \code{TRUE}, scale and offset values are ignored}
\item{crs}{character. Description of the Coordinate Reference System (map projection) in \code{PROJ.4}, \code{WKT} or \code{authority:code} notation. If this argument is missing, and the x coordinates are within -360 .. 360 and the y coordinates are within -90 .. 90, longitude/latitude is assigned}
\item{extent}{\code{\link{SpatExtent}}}
\item{...}{additional \code{SpatRaster} objects} 
}


\value{
SpatRasterDataset
}

\seealso{
\code{\link{describe}}
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

