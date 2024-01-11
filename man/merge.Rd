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
Merge multiple SpatRasters to create a new SpatRaster with a larger spatial extent. The SpatRasters must have the same origin and spatial resolution. In areas where the SpatRasters overlap, the values of the SpatRaster that is first in the sequence of arguments (or in the SpatRasterCollection) will be retained (unless \code{first=FALSE}). Cells that are \code{NA} are ignored unless \code{na.rm=FALSE}

There is also a method for merging SpatVector with a data.frame; that is, to join the data.frame to the attribute table of the SpatVector.
}

\usage{
\S4method{merge}{SpatRaster,SpatRaster}(x, y, ..., first=TRUE, na.rm=TRUE, filename="", overwrite=FALSE, wopt=list())

\S4method{merge}{SpatRasterCollection,missing}(x, first=TRUE, na.rm=TRUE, filename="", ...)

\S4method{merge}{SpatVector,data.frame}(x, y, ...)
}

\arguments{
  \item{x}{SpatRaster, SpatRasterCollection, or SpatVector}
  \item{y}{missing if \code{x} is a SpatRasterCollection. SpatRaster if \code{x} is a SpatRaster. data.frame if \code{x} is a SpatVector}
  \item{...}{if \code{x} is a SpatRaster: additional objects of the same class as \code{x}. If \code{x} is a SpatRasterCollection: options for writing files as in \code{\link{writeRaster}}. If \code{x} is a SpatVector, the same arguments as in \code{\link[base]{merge}}}
  \item{first}{logical. If \code{TRUE}, in areas where rasters overlap, the first value is used. Otherwise the last value is used}
  \item{na.rm}{logical. If \code{TRUE} values that are \code{NA} are ignored}
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

# if you have many SpatRasters, it may be convenient
# to make a SpatRasterCollection
rlist <- list(x, y, z)
rsrc <- sprc(rlist)

m <- merge(rsrc)


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

