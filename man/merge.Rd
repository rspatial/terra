\name{merge}

\docType{methods}

\alias{merge}
\alias{merge,SpatRaster,SpatRaster-method}
\alias{merge,SpatRasterCollection,missing-method}
\alias{merge,SpatExtent,SpatExtent-method}
\alias{merge,SpatVector,SpatVector-method}
\alias{merge,SpatVector,data.frame-method}

\title{ 
Merge multiple SpatRaster objects or SpatExtent objects, or merge a SpatVector with a data.frame
}

\description{
Merge multiple SpatRasters to create a new SpatRaster object with a larger spatial extent. The SpatRasters must have the same origin and spatial resolution. In areas where the SpatRasters overlap, the values of the SpatRaster that is first in the sequence of arguments (or in the SpatRasterCollection) will be retained (unless \code{first=FALSE}. 

See \code{\link{classify}} to merge a \code{SpatRaster} and a \code{data.frame}. You can also merge SpatExtent objects. 

There is a also a method for merging SpatVector with a data.frame; that is, to join the data.frame to the attribute table of the SpatVector.
}

\usage{
\S4method{merge}{SpatRaster,SpatRaster}(x, y, ..., first=TRUE, filename="", overwrite=FALSE, wopt=list())

\S4method{merge}{SpatRasterCollection,missing}(x, first=TRUE, filename="", ...)

\S4method{merge}{SpatExtent,SpatExtent}(x, y, ...)

\S4method{merge}{SpatVector,data.frame}(x, y, ...)
}

\arguments{
  \item{x}{SpatRaster or SpatExtent}
  \item{y}{object of same class as \code{x}}
  \item{...}{if \code{x} is a SpatRaster: additional objects of the same class as \code{x}. If \code{x} is a SpatRasterCollection: options for writing files as in \code{\link{writeRaster}}. If \code{x} is a SpatVector, the same arguments as in \code{\link[base]{merge}}}
  \item{first}{logical. If \code{TRUE}, in areas where rasters overlap, the first value that is not \code{NA} is used. Otherwise the last value that is not \code{NA} is used}
  \item{filename}{character. Output filename}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list with named options for writing files as in \code{\link{writeRaster}}}
} 


\value{
SpatRaster or SpatExtent
}

\seealso{
Combining tiles with \code{\link{vrt}} may be more efficient. See \code{\link{mosaic}} for averaging overlapping regions. 
}

\note{
You can use \code{merge} with \code{do.call} to merge a list of SpatRasters (see example). But note that if the list is named, these names are used by \code{merge}. So if all elements are named, there should be one element with a SpatRaster called \code{x} and another one called \code{y}. For example with \code{names(x)[1:2] <- c("x"m "y")}. You can also removed the names of the the first two elements (assuming these are SpatRasters) with \code{names(x)[1:2] <- ""}.
}


\examples{
x <- rast(xmin=-110, xmax=-80, ymin=40, ymax=70, ncols=30, nrows=30, vals=1)
y <- rast(xmin=-85, xmax=-55, ymax=60, ymin=30, ncols=30, nrows=30, vals=2)
z <- rast(xmin=-60, xmax=-30, ymax=50, ymin=20, ncols=30, nrows=30, vals=3)

m1 <- merge(x, y, z)
m2 <- merge(z, y, x)
m3 <- merge(y, x, z)

# if you have many SpatRasters make a SpatRasterCollection from a list
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

