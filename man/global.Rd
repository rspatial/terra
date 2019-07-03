\name{global}

\alias{global}

\alias{global,SpatRaster-method}

\title{global statistics}

\description{
Compute global statistics, that is summarized values of an entire SpatRaster. 

If \code{stat} is a true \code{function}, \code{global} will fail (gracefully) for very large SpatRasters, but it will in most cases work for functions that can be defined as by a character argument ("mean", "min", "max", or "sum"). 

If a function is used, it should accept a \code{na.rm} argument (or at least a \code{...} argument)
}

\usage{
\S4method{global}{SpatRaster}(x, fun="mean", na.rm=TRUE,  ...) 
}

\arguments{
  \item{x}{SpatRaster}
  \item{fun}{function to be applied to summarize the values by zone. Either as character: 'mean', 'sd', 'min', 'max', 'sum'; or, for relatively small SpatRasters, a proper function}
  \item{na.rm}{logical. If \code{TRUE}, \code{NA} values in \code{x} are ignored}
  \item{...}{Additional arguments. None implemented}  
}

\value{
A matrix with a value for each layer
}


\seealso{ See \code{\link{zonal}} for "zonal" statistics"  and \code{\link{reduce}} or \code{\link{local}} for "local" statistics, and \code{\link{extract}} for summarizing values for polygons}


\examples{
r <- rast(ncols=10, nrows=10)
values(r) <- 1:ncell(r)
global(r, "sum")
global(r, "mean")
}

\keyword{spatial}

