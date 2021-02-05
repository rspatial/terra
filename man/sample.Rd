\name{spatSample}

\docType{methods}

\alias{spatSample}
\alias{spatSample,SpatRaster-method}
\alias{spatSample,SpatExtent-method}
\alias{spatSample,SpatVector-method}

\title{Take a regular sample}

\description{
Take a spatial sample from a SpatRaster, SpatVector or SpatExtent. Sampling a SpatVector or SpatExtent always returns a SpatVector of points. 

With a SpatRaster, you can get cell values, cell numbers (\code{cells=TRUE}) or (when \code{type="regular"} and \code{as.raster=TRUE}) get a new SpatRaster with the same extent, but fewer cells. Note that, in order to assure regularity when requesting a regular sample, the number of values returned may not be exactly the same as the \code{size} requested.
}

\usage{
\S4method{spatSample}{SpatRaster}(x, size, method="regular", replace=FALSE, 
       na.rm=FALSE, as.raster=FALSE, cells=FALSE, ...)

\S4method{spatSample}{SpatVector}(x, size, method="regular", by_geom=TRUE, 
       strata=NULL, chess="", ...)

\S4method{spatSample}{SpatExtent}(x, size, method="regular", lonlat, ...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{size}{numeric. The sample size}
  \item{method}{character. Should be "regular" or "random". It can also be "stratified" if \code{x} is a SpatVector}
  \item{replace}{logical. If \code{TRUE}, sampling is with replacement (if \code{method="random"}}
  \item{na.rm}{logical. If \code{TRUE}, code{NAs} are removed. Only used with random sampling of cell values. That is with \code{method="random", as.raster=FALSE, cells=FALSE}}
  \item{as.raster}{logical. If \code{TRUE}, a SpatRaster is returned}
  \item{cells}{logical. If \code{TRUE}, cellnumbers are returned instead of values}
  \item{by_geom}{logical. If \code{TRUE}, sampling is done seperately for each geometry (multi-polygon or multi-line), otherwise, a single sample is taken for all geometries}
  \item{strata}{if not NULL, stratified random sampling is done, taking \code{size} samples from each stratum. If \code{x} has polygon geometry, \code{strata} must be a field name (or index) in \code{x}. If \code{x} has point geometry, \code{strata} can be a SpatVector of polygons or a SpatRaster}
  \item{chess}{character. One of "", "white", or "black". For stratified sampling if \code{strata} is a SpatRaster. If not "", samples are only taken from alternate cells, organized like the "white" or "black" fields on a chessboard}
  
  \item{lonlat}{logical. If \code{TRUE}, sampling of a SpatExtent is weighted by \code{cos(latitude)}. For SpatRaster and SpatVector this done based on the \code{\link{crs}}, but it is ignored if \code{as.raster=TRUE}}
  \item{...}{additional arguments. None implemented}
}


\details{
In some cases you may want to know the location of the sampled cells. In that situation you can take a sample of the cell numbers and use extract. See examples.

In stead of \code{spatSample(x, size, method="random")} you can also use the equivalent base method \code{sample(x, size)}. The base method also works for sampling the geometries of \code{SpatVector} (you can take a sample from the number of geometries and use that as an index). 
}

\value{
numeric or SpatRaster
}

\examples{
f <- system.file("ex/elev.tif", package="terra")
r <- rast(f)
s <- spatSample(r, 10, as.raster=TRUE)
spatSample(r, 10)
spatSample(r, 10, "random")

## if you require cell numbers and/or coordinates
size <- 6
# random cells 
cells <- spatSample(r, 6, "random", cells=TRUE)
v <- r[cells]
xy <- xyFromCell(r, cells)
cbind(xy, v)

# regular 
cells <- spatSample(r, 6, "regular", cells=TRUE)
v <- r[cells]
xy <- xyFromCell(r, cells)
cbind(xy, v)

## SpatExtent 
e <- ext(r)
spatSample(e, 10, "random", lonlat=TRUE)


## SpatVector
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
#sample geometries 
i <- sample(nrow(v), 5)
vv <- v[i,]

}

\keyword{spatial}
