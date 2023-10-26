\name{rasterizeGeom}

\docType{methods}

\alias{rasterizeGeom}
\alias{rasterizeGeom,SpatVector,SpatRaster-method}

\title{Rasterize geometric properties of vector data}

\description{
Rasterization of geometric properties of vector data. You can get the count of the number of geometries in each cell; the area covered by polygons; the length of the lines; or the number of lines that cross the boundary of each cell. See \code{\link{rasterize}} for standard rasterization (of attribute values associated with geometries).

The area of polygons is intended for summing the area of polygons that are relatively small relative to the raster cells, and for when there may be multiple polygons per cell. See \code{rasterize(fun="sum")} for counting large polygons and \code{rasterize(cover=TRUE)} to get the fraction that is covered by larger polygons.
}

\usage{
\S4method{rasterizeGeom}{SpatVector,SpatRaster}(x, y, fun="count", unit="m", filename="", ...)
}

\arguments{
  \item{x}{SpatVector}
  \item{y}{SpatRaster}
  \item{fun}{character. "count", "area", "length", or "crosses"}
  \item{unit}{character. "m" or "km"}
  \item{filename}{character. Output filename}
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
}

\seealso{
\code{\link{rasterize}}
}


\value{
SpatRaster
}

\examples{
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
r <- rast(v, res=.1)

# length of lines
lns <- as.lines(v)
x <- rasterizeGeom(lns, r, fun="length", "km")

# count of points
set.seed(44)
pts <- spatSample(v, 100)
y <- rasterizeGeom(pts, r)

# area of polygons
pols <- buffer(pts, 1000)
z <- rasterizeGeom(pols, r, fun="area")
}

\keyword{spatial}

