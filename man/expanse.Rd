\name{expanse}

\alias{expanse}
\alias{expanse,SpatRaster-method}
\alias{expanse,SpatVector-method}

  
\title{Get the expanse (area) of individual polygons or for all (summed) raster cells}

\description{
Compute the area covered by polygons or for all raster cells that are not \code{NA}. 

This method computes areas for longitude/latitude rasters, as the size of the cells is constant in degrees, but not in square meters. But it can also be important if the coordinate reference system is planar, but not equal-area. 

For vector data, the best way to compute area is to use the longitude/latitude CRS. This is contrary to (erroneous) popular belief that suggest that you should use a planar coordinate reference system. This is done automatically, if \code{transform=TRUE}.
}

\usage{
\S4method{expanse}{SpatRaster}(x, unit="m", transform=TRUE, byValue=FALSE)

\S4method{expanse}{SpatVector}(x, unit="m", transform=TRUE)
}

\arguments{
  \item{x}{SpatRaster or SpatVector}
  \item{unit}{character. One of "m", "km", or "ha"}
  \item{transform}{logical. If \code{TRUE}, planar CRS are transformed to lon/lat for accuracy}
  \item{byValue}{logical. If \code{TRUE}, the area for each unique cell value is returned}
}

\value{
numeric.  If \code{x} has no values, the total size of all cells. Otherwise, the total area size of all cells that are not \code{NA}, expressed in square meters, square kilometers, or hectares.

If \code{byValue=TRUE} a matrix is returned with three columns (layer, value, area) 
}


\seealso{
\code{\link{cellSize}} for a the size of individual cells of a raster, that can be summed with \code{\link{global}} or with \code{\link{zonal}} to get the area for different categories.
}

\examples{

### SpatRaster 
r <- rast(nrows=18, ncols=36)
v <- 1:ncell(r)
v[200:400] <- NA
values(r) <- v

# summed area in km2
expanse(r, unit="km")

# all cells 
expanse(rast(r), unit="km")

r <- rast(ncols=90, nrows=45, ymin=-80, ymax=80)
m <- project(r, "+proj=merc")

expanse(m, unit="km")
expanse(m, unit="km", transform=FALSE)

m2 <- c(m, m)
values(m2) <- cbind(c(1,2,NA,NA), c(11:14))
expanse(m2, unit="km", byValue=TRUE)



### SpatVector
v <- vect(system.file("ex/lux.shp", package="terra"))

a <- expanse(v)
a
sum(a)
}

\keyword{methods}
\keyword{spatial}
