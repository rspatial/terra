\name{disagg}

\docType{methods}

\alias{disagg}
\alias{disagg,SpatRaster-method}
\alias{disagg,SpatVector-method}

\title{Disaggregate raster cells or vector geometries}

\description{
\code{SpatRaster}: Create a SpatRaster with a higher resolution (smaller cells). The values in the new SpatRaster are the same as in the larger original cells.

\code{SpatVector}: Separate multi-objects (points, lines, polygons) into single objects; or further into segments (for lines or polygons).
}

\usage{
\S4method{disagg}{SpatRaster}(x, fact, method="near", filename="", ...)

\S4method{disagg}{SpatVector}(x, segments=FALSE)
}

\arguments{
  \item{x}{SpatRaster or SpatVector}
  \item{fact}{positive integer. Aggregation factor expressed as number of cells in each direction (horizontally and vertically). Or two integers (horizontal and vertical aggregation factor) or three integers (when also aggregating over layers)}  
  \item{method}{character. Either "near" for nearest or "bilinear" for bilinear interpolation}
  \item{segments}{logical. Should (poly-)lines or polygons be disaggregated into their line-segments?}
  \item{filename}{character. Output filename}
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
}


\seealso{\code{\link{aggregate}}, \code{\link{resample}}}

\value{
SpatRaster
}


\examples{
r <- rast(ncols=10, nrows=10)
rd <- disagg(r, fact=c(10, 2))
ncol(rd)
nrow(rd)
values(r) <- 1:ncell(r)
rd <- disagg(r, fact=c(4, 2))

}


\keyword{methods}
\keyword{spatial}

