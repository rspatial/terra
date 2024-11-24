\name{surfArea}

\alias{surfArea}
\alias{surfArea,SpatRaster-method}

\title{
Compute surface area from elevation data
}

\description{
It is often said that if Wales was flattened out it would have an area bigger than England. This function computes the surface area for a raster with elevation values, taking into account the sloping nature of the surface.
}

\usage{
surfArea(x, filename="", ...)
}

\arguments{
  \item{x}{SpatRaster with elevation values. Currently the raster CRS must be planar and have the same distance units (e.g. m) as the elevation values}
  \item{filename}{character. Output filename}
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
}

\value{SpatRaster}

\references{
Jenness, Jeff S., 2004. Calculating Landscape Surface Area from Digital Elevation Models. Wildlife Society Bulletin 32(3): 829-839
}

\author{
Barry Rowlingson <b.rowlingson@lancaster.ac.uk>
}

\examples{
v <- rast(volcano, crs="local")
x <- terra::surfArea(v)
}

\keyword{spatial}

