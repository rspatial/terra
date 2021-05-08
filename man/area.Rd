\name{size}

\alias{size}
\alias{size,SpatRaster-method}
\alias{size,SpatVector-method}

\alias{area}
\alias{area,SpatRaster-method}
\alias{area,SpatVector-method}

  
\title{Get the size (area) for individual polygons or for all (summed) raster cells}

\description{
Compute the area covered by polygons or for all raster cells that are not \code{NA}. 

This method computes areas for longitude/latitude rasters, as the size of the cells is constant in degrees, but not in meters. But it can also be important if the coordinate reference system is planar, but not equal-area. 

For vector data, the best way to compute area is to use the longitude/latitude CRS. This is contrary to (erroneous) popular belief that suggest that you should use a planar coordinate reference system. This is done automatically, if \code{transform=TRUE}.
}

\usage{
\S4method{size}{SpatRaster}(x, unit="m", transform=TRUE)

\S4method{size}{SpatVector}(x, unit="m", transform=TRUE)
}


\arguments{
  \item{x}{SpatRaster or SpatVector}
  \item{unit}{character. One of "m", "km", or "ha"}
  \item{transform}{logical. If \code{TRUE}, planar CRS are transformed to lon/lat for accuracy}
}

\value{
numeric 
}

\seealso{\code{\link{cellSize}} for a the size of individual cells of a raster}

\examples{

### SpatRaster 
r <- rast(nrow=18, ncol=36)
v <- 1:ncell(r)
v[200:400] <- NA
values(r) <- v

# summed area in km2
size(r, unit="km")


r <- rast(ncol=90, nrow=45, ymin=-80, ymax=80)
m <- project(r, "+proj=merc")

a <- suze(m, unit="km")
b <- global(!is.na(r), "sum")


### SpatVector
v <- vect(system.file("ex/lux.shp", package="terra"))

a <- size(v)
a
sum(a)
}

\keyword{methods}
\keyword{spatial}
