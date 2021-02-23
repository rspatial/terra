\name{area}

\alias{area}
\alias{perimeter}
\alias{area,SpatRaster-method}
\alias{area,SpatVector-method}
\alias{perimeter,SpatVector-method}
  
\title{Area and perimeter}

\description{
Compute the area of polygons or for raster cells that are not \code{NA}. Computing the surface area of raster cells is particularly relevant for longitude/latitude rasters, as the size of the cells is constant in degrees, but not in meters. But it can also be very important with raster data if the coordinate reference system us not equal-area. In that case, you can use the \code{correct=TRUE} option. 

For vector data, the best way to compute area is to use the longitude/latitude crs if that is what the data come in. This is contrary to (erroneous) popular belief that suggest that you should use a planar crs.

The perimeter method works only on SpatVector objects, and computes the length of lines or the perimeter of polygons.
}

\usage{
\S4method{area}{SpatRaster}(x, sum=TRUE, correct=FALSE, filename="", overwrite=FALSE, ...)

\S4method{area}{SpatVector}(x, ...)

\S4method{perimeter}{SpatVector}(x)
}


\arguments{
  \item{x}{SpatRaster or SpatVector}
  \item{sum}{logical. If \code{TRUE} the summed area of the cells that are not \code{NA} is returned. Otherwise, a SpatRaster with the area for each cell is returned}
  \item{correct}{logical. If \code{TRUE}, the area is not computed based on the linear units of the coordinate reference system, but on the *actual* area, correcting for distortion. This may be considerably slower}
  \item{filename}{character. Output filename. Optional}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
}

\value{
area or perimeter in m2
}

\examples{

### SpatRaster 
r <- rast(nrow=18, ncol=36)
v <- 1:ncell(r)
v[200:400] <- NA
values(r) <- v

# area for each raster cell
a <- area(r, sum=FALSE)

# summed area in km2
area(r) / 1000000

## you can use mask to remove the cells in r that are NA
## and compute the global sum to get the same result
#am <- mask(a, r)
#global(am, "sum", na.rm=TRUE) / 1000000

# effect of "correct"
#r <- rast(ncol=90, nrow=45, ymin=-80, ymax=80)
#m <- project(r, "+proj=merc")

#a <- area(m, sum=FALSE, wopt=list(names="naive"))
#b <- area(m, correct=TRUE, sum=FALSE, names="corrected")
#plot(c(a, b)/1000000, nc=1)


### SpatVector
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
a <- area(v)
a
sum(a)
perimeter(v)

}

\keyword{methods}
\keyword{spatial}
