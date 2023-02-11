\name{viewshed}

\alias{viewshed}
\alias{viewshed,SpatRaster-method}

\title{Compute a viewshed}

\description{
Use elevation data to compute the locations that can be seen, or how much higher they would have to be to be seen, from a certain position. The raster data coordinate reference system must planar (not lon/lat), with the elevation values in the same unit as the distance unit of the coordinate reference system. 
}

\usage{
\S4method{viewshed}{SpatRaster}(x, loc, observer=1.80, target=0, curvcoef=6/7, output="yes/no", filename="", ...) 
}

\arguments{
  \item{x}{SpatRaster, single layer with elevation values. Values should have the same unit as the map units}
  \item{loc}{location (x and y coordinates) or a cell number}
  \item{observer}{numeric. The height above the elevation data of the observer}
  \item{target}{numeric. The height above the elevation data of the targets}
  \item{curvcoef}{numeric. Coefficient to consider the effect of the curvature of the earth and refraction of the atmosphere. The elevation values are corrected with: \code{elevation = elevation - curvcoeff * (distance)^2 / (earth_diameter)}. This means that with the default value of 0.85714, you lose sight of about 1 meter of elevation for each 385 m of planar distance}
  \item{output}{character. Can be "yes/no" to get a binary (logical) output showing what areas are visible; "land" to get the height above the current elevation that would be visible; or "sea" the elevation above sea level that would be visible}
  \item{filename}{character. Output filename}
  \item{...}{Options for writing files as in \code{\link{writeRaster}}}
}

\seealso{\code{\link{terrain}}}

\references{
The algorithm used is by Wang et al.: https://www.asprs.org/wp-content/uploads/pers/2000journal/january/2000_jan_87-90.pdf.
}

\examples{

if (gdal() >= 3.1) {
	f <- system.file("ex/elev.tif", package="terra")
	r <- rast(f)
	x <- project(r, "EPSG:2169")
	p <- cbind(70300, 96982)
	v <- viewshed(x, p, 0, 0, 0.85714)
}

}

\keyword{spatial}
