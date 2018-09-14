\name{aggregate}

\docType{methods}

\alias{aggregate}
\alias{aggregate,SpatRaster-method}

\title{Aggregate raster cells}

\description{
Aggregate a SpatRaster object to create a new RasterLayer or RasterBrick with a lower resolution (larger cells). 
Aggregation groups rectangular areas to create larger cells. The value for the resulting cells is computed with a user-specified function. 
}

\usage{
\S4method{aggregate}{SpatRaster}(x, fact=2, fun="mean", na.rm=TRUE, filename='', ...)
}

\arguments{
  \item{x}{Raster* object or SpatialPolygons* object}
  \item{fact}{integer. Aggregation factor expressed as number of cells in each direction (horizontally and vertically). Or two integers (horizontal and vertical aggregation factor) or three integers (when also aggregating over layers)}  
  \item{fun}{function used to aggregate values }  
  \item{na.rm}{logical. If \code{TRUE}, NA cells are removed from calculations }
  \item{filename}{character. Output filename (optional)}
  \item{...}{if \code{x} is a Raster* object, additional arguments as for \code{\link{writeRaster}}}
}


\details{
Aggregation starts at the upper-left end of a raster. If a division of the number of columns or rows with \code{factor} does not return an integer, the extent of the resulting Raster object will be somewhat larger then that of the original SpatRaster. For example, if an input RasterLayer has 100 columns, and \code{fact=12}, the output Raster object will have 9 columns and the maximum x coordinate of the output SpatRaster is also adjusted.
  
The function \code{fun} should take multiple numbers, and return a single number. For example \code{mean}, \code{modal}, \code{min} or \code{max}. 
It should also accept a \code{na.rm} argument (or ignore it as one of the 'dots' arguments). 
}

\value{
SpatRaster object
}


\examples{
r <- georst()
# a new aggregated raster, no values
ra <- aggregate(r, fact=10)

values(r) <- runif(ncell(r))

# a new aggregated raster, max of the values
ra <- aggregate(r, fact=10, fun=max)

\dontrun{
# multiple layers
s <- stack(r, r*2)
x <- aggregate(s,2)
}
}

\keyword{methods}
\keyword{spatial}

