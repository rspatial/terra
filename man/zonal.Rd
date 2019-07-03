\name{zonal}

\alias{zonal}

\alias{zonal,SpatRaster,SpatRaster-method}

\title{Zonal statistics}

\description{
Compute zonal statistics, that is summarized values of a SpatRaster for each "zone" defined by another SpatRaster. 

If \code{stat} is a true \code{function}, \code{zonal} will fail (gracefully) for very large Raster objects, but it will in most cases work for functions that can be defined as by a character argument ('mean', 'sd', 'min', 'max', or 'sum'). In addition you can use 'count' to count the number of cells in each zone (only useful with \code{na.rm=TRUE}, otherwise \code{freq(z)} would be more direct.

If a function is used, it should accept a \code{na.rm} argument (or at least a \code{...} argument)
}

\usage{
\S4method{zonal}{SpatRaster,SpatRaster}(x, z, fun="mean", na.rm=TRUE, ...) 
}

\arguments{
  \item{x}{SpatRaster}
  \item{z}{SpatRaster with values representing zones}
  \item{fun}{function to be applied to summarize the values by zone. Either as character: 'mean', 'sd', 'min', 'max', 'sum'; or, for relatively small SpatRasters, a proper function}
  \item{na.rm}{logical. If \code{TRUE}, \code{NA} values in \code{x} are ignored}
  \item{...}{Additional arguments. None implemented}  
}

\value{
A matrix with a value for each zone (unique value in \code{zones})
}


\seealso{ See \code{\link{global}} for "global" statistics (i.e., all of \code{x} is considered a single zone),  \code{\link{reduce}} for local statistics, and \code{\link{extract}} for summarizing values for polygons}



\examples{
r <- rast(ncols=10, nrows=10)
values(r) <- 1:ncell(r)
z <- rast(r)
values(z) <- rep(1:4, each=25)
zonal(r, z, "sum")
}

\keyword{spatial}

