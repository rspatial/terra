\name{tapp}

\docType{methods}

\alias{tapp}
\alias{tapp,SpatRaster-method}

\title{Apply a function to subsets of layers of a SpatRaster}

\description{
Apply a function to subsets of layers of a SpatRaster (similar to \code{\link[base]{tapply}} and \code{\link[stats]{aggregate}}). The layers are combined based on the \code{index}.

The number of layers in the output SpatRaster equals the number of unique values in \code{index} times the number of values that the supplied function returns for a single vector of numbers.

For example, if you have a SpatRaster with 6 layers, you can use \code{index=c(1,1,1,2,2,2)} and \code{fun=sum}. This will return a SpatRaster with two layers. The first layer is the sum of the first three layers in the input SpatRaster, and the second layer is the sum of the last three layers in the input SpatRaster. Indices are recycled such that \code{index=c(1,2)} would also return a SpatRaster with two layers (one based on the odd layers (1,3,5), the other based on the even layers (2,4,6)). 

The index can also be one of the following values to group by time period (if \code{x} has the appropriate \code{\link{time}} values): "years", "months", "yearmonths", "week" (the ISO 8601 week number, see Details), "yearweeks", "days", "doy" (day of the year), "7days" (seven-day periods starting at Jan 1 of each year), "10days", or "15days". It can also be a function that makes groups from time values.

See \code{\link{app}} or \code{\link{Summary-methods}} if you want to use a more efficient function that returns multiple layers based on \bold{all} layers in the SpatRaster.
}

\usage{
\S4method{tapp}{SpatRaster}(x, index, fun, ..., cores=1, filename="", overwrite=FALSE, wopt=list())
}

\arguments{
  \item{x}{SpatRaster}
  \item{index}{factor or numeric (integer). Vector of length \code{nlyr(x)} (shorter vectors are recycled) grouping the input layers. It can also be one of the following values: "years", "months", "yearmonths", "days", "week" (ISO 8601 week number), or "doy" (day of the year)}
  \item{fun}{function to be applied. The following functions have been re-implemented in C++ for speed: "sum", "mean", "median", "modal", "which", "which.min", "which.max", "min", "max", "prod", "any", "all", "sd", "std", "first". To use the base-R function for say, "min", you could use something like \code{fun = \(i) min(i)}}
  \item{...}{additional arguments passed to \code{fun}}
  \item{cores}{positive integer. If \code{cores > 1}, a 'parallel' package cluster with that many cores is created and used. You can also supply a cluster object. Ignored for functions that are implemented by terra in C++ (see under fun)}  
  \item{filename}{character. Output filename}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list with named options for writing files as in \code{\link{writeRaster}}}
}

\value{
SpatRaster
}

\details{
"week" follows the ISO 8601 definition. Weeks start on Monday. If the week containing 1 January has four or more days in the new year, then it is considered week "01". Otherwise, it is the last week of the previous year (week "52" or "53", and the next week is week 1.
}

\seealso{\code{\link{app}}, \code{\link{Summary-methods}}}

\examples{
r <- rast(ncols=10, nrows=10)
values(r) <- 1:ncell(r)
s <- c(r, r, r, r, r, r)
s <- s * 1:6
b1 <- tapp(s, index=c(1,1,1,2,2,2), fun=sum)
b1
b2 <- tapp(s, c(1,2,3,1,2,3), fun=sum)
b2
}

\keyword{methods}
\keyword{spatial}

