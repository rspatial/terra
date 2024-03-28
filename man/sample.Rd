\name{spatSample}

\docType{methods}

\alias{spatSample}
\alias{spatSample,SpatRaster-method}
\alias{spatSample,SpatExtent-method}
\alias{spatSample,SpatVector-method}

\title{Take a regular sample}

\description{
Take a spatial sample from a SpatRaster, SpatVector or SpatExtent. Sampling a SpatVector or SpatExtent always returns a SpatVector of points. 

With a SpatRaster, you can get cell values, cell numbers (\code{cells=TRUE}), coordinates (\code{xy=TRUE}) or (when \code{method="regular"} and \code{as.raster=TRUE}) get a new SpatRaster with the same extent, but fewer cells. 

In order to assure regularity when requesting a regular sample, the number of cells or points returned may not be exactly the same as the \code{size} requested.
}

\usage{
\S4method{spatSample}{SpatRaster}(x, size, method="random", replace=FALSE, na.rm=FALSE, 
    as.raster=FALSE, as.df=TRUE, as.points=FALSE, values=TRUE, cells=FALSE, 
    xy=FALSE, ext=NULL, warn=TRUE, weights=NULL, exp=5, exhaustive=FALSE)

\S4method{spatSample}{SpatVector}(x, size, method="random", strata=NULL, chess="")

\S4method{spatSample}{SpatExtent}(x, size, method="random", lonlat, as.points=FALSE)
}

\arguments{
  \item{x}{SpatRaster, SpatVector or SpatExtent}
  \item{size}{numeric. The sample size. If \code{x} is a SpatVector, you can also provide a vector of the same length as \code{x} in which case sampling is done separately for each geometry. If \code{x} is a SpatRaster, and you are using \code{method="regular"} you can specify the size as two numbers (number of rows and columns)}
  \item{method}{character. Should be "regular" or "random", If \code{x} is a \code{SpatRaster}, it can also be "stratified" (each value in \code{x} is a stratum) or "weights" (each value in \code{x} is a probability weight)}
  \item{replace}{logical. If \code{TRUE}, sampling is with replacement (if \code{method="random"})}
  \item{na.rm}{logical. If \code{TRUE}, \code{NAs} are removed. Only used with random sampling of cell values. That is with \code{method="random", as.raster=FALSE, cells=FALSE}}
  \item{as.raster}{logical. If \code{TRUE}, a SpatRaster is returned}
  \item{as.df}{logical. If \code{TRUE}, a data.frame is returned instead of a matrix}
  \item{as.points}{logical. If \code{TRUE}, a SpatVector of points is returned}
  \item{values}{logical. If \code{TRUE} raster cell values are returned}
  \item{cells}{logical. If \code{TRUE}, cell numbers are returned. If \code{method="stratified"} this is always set to \code{TRUE} if \code{xy=FALSE}}
  \item{xy}{logical. If \code{TRUE}, cell coordinates are returned}
  \item{ext}{SpatExtent or NULL to restrict sampling to a subset of the area of \code{x}}
  \item{warn}{logical. Give a warning if the sample size returned is smaller than requested}
  \item{weights}{SpatRaster. Used to provide weights when \code{method="stratified"}}

  \item{strata}{if not NULL, stratified random sampling is done, taking \code{size} samples from each stratum. If \code{x} has polygon geometry, \code{strata} must be a field name (or index) in \code{x}. If \code{x} has point geometry, \code{strata} can be a SpatVector of polygons or a SpatRaster}
  \item{chess}{character. One of "", "white", or "black". For stratified sampling if \code{strata} is a SpatRaster. If not "", samples are only taken from alternate cells, organized like the "white" or "black" fields on a chessboard}
  
  \item{lonlat}{logical. If \code{TRUE}, sampling of a SpatExtent is weighted by \code{cos(latitude)}. For SpatRaster and SpatVector this done based on the \code{\link{crs}}, but it is ignored if \code{as.raster=TRUE}}
  
  \item{exp}{numeric >= 1. "Expansion factor" that is multiplied with \code{size} to get an initial sample used for stratified samples and random samples with \code{na.rm=TRUE} to try to get at least \code{size} samples}
  
  \item{exhaustive}{logical. If \code{TRUE} and \code{na.rm=TRUE} first all cells that are not \code{NA} are determined and a sample is taked from these cells. This is useful when you are dealing with a very large raster that is sparse (most cells are NA). Otherwise, the default approach may not find enough samples. This should not be used in other cases, especially not with large rasters that mostly have values}  
}



\value{
numeric matrix, data.frame, SpatRaster or SpatVector
}

\examples{
f <- system.file("ex/elev.tif", package="terra")
r <- rast(f)
s <- spatSample(r, 10, as.raster=TRUE)
spatSample(r, 5)
spatSample(r, 5, na.rm=TRUE)
spatSample(r, 5, "regular")

## if you require cell numbers and/or coordinates
size <- 6
spatSample(r, 6, "random", cells=TRUE, xy=TRUE, values=FALSE)

# regular, with values 
spatSample(r, 6, "regular", cells=TRUE, xy=TRUE)

# stratified
rr <- rast(ncol=10, nrow=10, names="stratum")
set.seed(1)
values(rr) <- round(runif(ncell(rr), 1, 3))
spatSample(rr, 2, "stratified", xy=TRUE)

s <- spatSample(rr, 5, "stratified", as.points=TRUE)
plot(rr, plg=list(title="raster"))
plot(s, 1, add=TRUE, plg=list(x=185, y=1, title="points"))
 
## SpatExtent 
e <- ext(r)
spatSample(e, 10, "random", lonlat=TRUE)

## SpatVector
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)

# sample the geometries 
i <- sample(v, 3)

# sample points in geometries
p <- spatSample(v, 3)

}

\keyword{spatial}
