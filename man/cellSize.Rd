\name{cellSize}

\alias{cellSize}
\alias{cellSize,SpatRaster-method}
  
\title{Area covered by each raster cell}

\description{

Compute the area covered by individual raster cells. 

Computing the surface area of raster cells is especially relevant for longitude/latitude rasters. 

But note that for both angular (longitude/latitude) and for planar (projected) coordinate reference systems raster cells sizes are generally not constant, unless you are using an equal-area coordinate reference system. 

For planar CRSs, the area is therefore not computed based on the linear units of the coordinate reference system, but on the *actual* area by transforming cells to longitude/latitude. If you do not want that correction, you can use \code{transform=FALSE} or \code{init(x, prod(res(x)))}
}

\usage{
\S4method{cellSize}{SpatRaster}(x, mask=FALSE, lyrs=FALSE, unit="m", transform=TRUE, rcx=100, filename="", ...)
}


\arguments{
  \item{x}{SpatRaster}
  \item{mask}{logical. If \code{TRUE}, cells that are \code{NA} in \code{x} are also \code{NA} in the output}
  \item{lyrs}{logical. If \code{TRUE} and \code{mask=TRUE}, the output has the same number of layers as \code{x}. That is only useful if cases where the layers of \code{x} have different cells that are \code{NA}}
  \item{unit}{character. One of "m", "km", or "ha"}
  \item{transform}{logical. If \code{TRUE}, planar CRS data are transformed to lon/lat for accuracy}
  \item{rcx}{positive integer. The maximum number of rows and columns to be used to compute area of planar data if \code{transform=TRUE}. If \code{x} has more rows and/or columns, the raster is aggregated to match this limit, and values for the original cells are estimated by bilinear interpolation (see \code{resample}). This can save a lot of time}
  \item{filename}{character. Output filename}
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
}

\value{
SpatRaster
}

\seealso{\code{\link{expanse}}}

\examples{

# SpatRaster 
r <- rast(nrows=18, ncols=36)
v <- 1:ncell(r)
v[200:400] <- NA
values(r) <- v

# size of each raster cell
a <- cellSize(r)

# illustration of distortion
r <- rast(ncols=90, nrows=45, ymin=-80, ymax=80)
m <- project(r, "+proj=merc")

bad <- init(m, prod(res(m)) / 1000000, names="naive")
good <- cellSize(m, unit="km", names="corrected")
plot(c(good, bad), nc=1, mar=c(2,2,1,6))
}

\keyword{methods}
\keyword{spatial}
