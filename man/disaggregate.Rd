\name{disaggregate}

\docType{methods}

\alias{disaggregate}
\alias{disaggregate,SpatRaster-method}

\title{Disaggregate raster cells}

\description{
Disaggregate a SpatRaster to create a new SpatRaster with a higher resolution (smaller cells). The values in the new SpatRaster are the same as in the larger original cells.
}

\usage{
\S4method{disaggregate}{SpatRaster}(x, fact, filename="", overwrite=FALSE, wopt=list(), ...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{fact}{positive integer. Aggregation factor expressed as number of cells in each direction (horizontally and vertically). Or two integers (horizontal and vertical aggregation factor) or three integers (when also aggregating over layers)}  
  \item{filename}{character. Output filename. Optional}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list. Options for writing files as in \code{\link{writeRaster}}}
  \item{...}{additional arguments. None implemented}
}


\seealso{\code{\link{aggregate}}}

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

