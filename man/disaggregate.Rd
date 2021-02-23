\name{disaggregate}

\docType{methods}

\alias{disaggregate}
\alias{disaggregate,SpatRaster-method}
\alias{disaggregate,SpatVector-method}

\title{Disaggregate raster cells}

\description{
\code{SpatRaster}: Create a SpatRaster with a higher resolution (smaller cells). The values in the new SpatRaster are the same as in the larger original cells.

\code{SpatVector}: Separate multi-objects (points, lines, polygons) into single objects.
}

\usage{
\S4method{disaggregate}{SpatRaster}(x, fact, method="near", filename="", ...)

\S4method{disaggregate}{SpatVector}(x)
}

\arguments{
  \item{x}{SpatRaster or SpatVector}
  \item{fact}{positive integer. Aggregation factor expressed as number of cells in each direction (horizontally and vertically). Or two integers (horizontal and vertical aggregation factor) or three integers (when also aggregating over layers)}  
  \item{method}{character. Either "near" for nearest or "bilinear" for bilinear interpolation}
  \item{filename}{character. Output filename}
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
}


\seealso{\code{\link{aggregate}}, \code{\link{resample}}}

\value{
SpatRaster
}


\examples{
r <- rast(ncols=10, nrows=10)
rd <- disaggregate(r, fact=c(10, 2))
ncol(rd)
nrow(rd)
values(r) <- 1:ncell(r)
rd <- disaggregate(r, fact=c(4, 2))

}


\keyword{methods}
\keyword{spatial}

