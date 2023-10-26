\name{global}

\alias{global}

\alias{global,SpatRaster-method}

\title{global statistics}

\description{
Compute global statistics, that is summarized values of an entire SpatRaster. 

If \code{x} is very large \code{global} can fail, except when \code{fun} is one of these built-in functions "mean", "min", "max", "sum", "prod", "range" (min and max), "rms" (root mean square), "sd" (sample standard deviation), "std" (population standard deviation), "isNA" (number of cells that are NA), "notNA" (number of cells that are not NA).

The reason that this can fail with large raster and a custom function is that all values need to be loaded into memory. To circumvent this problem you can run \code{global} with a sample of the cells.

You can compute a weighted mean or sum by providing a SpatRaster with weights.
}

\usage{
\S4method{global}{SpatRaster}(x, fun="mean", weights=NULL, maxcell=Inf, ...) 
}

\arguments{
  \item{x}{SpatRaster}
  \item{fun}{function to be applied to summarize the values by zone. Either as one or more of these built-in character values: "max", "min", "mean", "sum", "range", "rms" (root mean square), "sd", "std" (population sd, using \code{n} rather than \code{n-1}), "isNA", "notNA"; or a proper R function (but these may fail for very large SpatRasters unless you specify \code{maxcell})}
  \item{...}{additional arguments passed on to \code{fun}}  
  \item{weights}{NULL or SpatRaster}  
  \item{maxcell}{positive integer used to take a regular sample of \code{x}. Ignored by the built-in functions.}  
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
x <- c(r, r/10)
global(x, c("sum", "mean", "sd"), na.rm=TRUE)

global(x, function(i) min(i) / max(i))
}

\keyword{spatial}

