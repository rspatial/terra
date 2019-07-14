\name{area}

\alias{area}
\alias{area,SpatRaster-method}
\alias{area,SpatVector-method}
\alias{length,SpatVector-method}
  
\title{Area and length}

\description{
Compute area of polygons or and the surfacer area or raster cells.
Computing the surface area of raster cells is only relevant for longitude/latitude rasters. For planar data, the cell area is constant, can be computed with \code{prod(res(x))}. For polygons, the area of each geometry is returned.

The length method works only on SpatVector objects, and computes the length of lines or the circumference of polygons.
}

\usage{
\S4method{area}{SpatRaster}(x, filename="", overwrite=FALSE,  wopt=list(), ...)

\S4method{area}{SpatVector}(x, ...)
\S4method{length}{SpatVector}(x)
}


\arguments{
  \item{x}{SpatRaster or SpatVector}
  \item{filename}{character. Output filename. Optional}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list. Options for writing files as in \code{\link{writeRaster}}}
  \item{...}{additional arguments. None implemented} 
}

\value{
If the coordinate reference system is longitude/latitude the values returned are in square meter for area and meter for length.

In other cases, the unit is set by coordinate references system. In most cases it is meter too. 
}



\examples{
f <- system.file("exdata/lux.shp", package="terra")
v <- vect(f)
a <- area(v)
a
sum(a)
length(v)
}

\keyword{methods}
\keyword{spatial}
