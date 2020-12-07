\name{rapp}

\docType{methods}

\alias{rapp}
\alias{rapp,SpatRaster-method}

\title{Apply a function to a range of the layers of a SpatRaster}

\description{
Apply a function to a range of the layers of a SpatRaster. The range is specified for each cell seperately by a two-layer SpatRaster \code{index}. 

The function used should return a single value.

See \code{\link{app}} or \code{\link{Summary-methods}} if you want to apply a function to all layers (or a subset of all layers) in a SpatRaster.
}

\usage{
\S4method{rapp}{SpatRaster}(x, index, fun, ..., filename="", overwrite=FALSE, wopt=list())
}

\arguments{
  \item{x}{SpatRaster}
  \item{index}{factor or numeric (integer). Vector of length \code{nlyr(x)} (shorter vectors are recycled) grouping the input layers}
  \item{fun}{function to be applied}
  \item{...}{additional arguments passed to \code{fun}}  
  \item{filename}{character. Output filename. Optional}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list. Options for writing files as in \code{\link{writeRaster}}}
}

\value{
SpatRaster
}

\seealso{\code{\link{app}}, \code{\link{Summary-methods}}, \code{\link{lapp}}, \code{\link{tapp}}}

\examples{
r <- rast(ncol=9, nrow=9)
values(r) <- 1:ncell(r)
s <- c(r, r, r, r, r, r)
s <- s * 1:6

start <- end <- rast(r)
start[] <- 1:3
end[]   <- 4:6
index   <- c(start, end)

rapp(s, index, fun="mean")
}

\keyword{methods}
\keyword{spatial}

