\name{viewshed}

\alias{viewshed}
\alias{viewshed,SpatRaster-method}

\title{Compute a viewshed}

\description{
Compute a locations that can be seen from a certain position (viewshed) based on elevation data. The raster data must be planar, with the elevation values in the same unit as the distance unit of the coordinate reference system. 
}

\usage{
\S4method{viewshed}{SpatRaster}(x, loc, observer=1.80, target=0, curvcoef=0.85714, filename="", ...) 
}

\arguments{
  \item{x}{SpatRaster, single layer with elevation values. Values should have the same unit as the map units}
  \item{loc}{location (x and y coordinates) or a cell number}
  \item{observer}{numeric. The height above the elevation data of the observer}
  \item{target}{numeric. The height above the elevation data of the targets}
  \item{curvcoef}{numeric. coefficient to correct for the curvature of the earth}
  \item{filename}{character. Output filename}
  \item{...}{Options for writing files as in \code{\link{writeRaster}}}
}

\seealso{\code{\link{viewshed}}}

\references{
The algorithm used is by Wang et al.: https://www.asprs.org/wp-content/uploads/pers/2000journal/january/2000_jan_87-90.pdf.
}

\examples{

f <- system.file("ex/elev.tif", package="terra")
r <- rast(f)
x <- project(r, "EPSG:2169")
v <- viewshed(x, c(70300, 96982), 0, 0, 0.85714)

}

\keyword{spatial}
