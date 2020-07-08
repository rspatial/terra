\name{global}

\alias{global}

\alias{global,SpatRaster-method}

\title{global statistics}

\description{
Compute global statistics, that is summarized values of an entire SpatRaster. 

If \code{x} is very large \code{global} will fail, except when \code{fun} is one of "mean", "min", "max", or "sum". 

You can compute a weighted mean or sum by providing a SpatRaster with weights.
}

\usage{
\S4method{global}{SpatRaster}(x, fun="mean", weights=NULL, ...) 
}

\arguments{
  \item{x}{SpatRaster}
  \item{fun}{function to be applied to summarize the values by zone. Either as character:"mean", "min", "max", "sum"; or, for relatively small SpatRasters, a proper function}
  \item{...}{additional arguments passed on to \code{fun}}  
  \item{weights}{NULL or SpatRaster}  
}

\value{
A \code{data.frame} with a row for each layer
}


\seealso{\code{\link{zonal}} for "zonal" statistics, and \code{\link{app}} or \code{\link{Summary-methods}} for "local" statistics, and \code{\link{extract}} for summarizing values for polygons. Also see \code{\link{focal}} for "focal" or "moving window" operations.}


\examples{
r <- rast(ncols=10, nrows=10)
values(r) <- 1:ncell(r)
global(r, "sum")
global(r, "mean", na.rm=TRUE)
}

\keyword{spatial}

