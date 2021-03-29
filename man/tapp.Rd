\name{tapp}

\docType{methods}

\alias{tapp}
\alias{tapp,SpatRaster-method}

\title{Apply a function to subsets of layers of a SpatRaster}

\description{
Apply a function to subsets of layers of a SpatRaster (similar to \code{\link[base]{tapply}} and \code{\link[stats]{aggregate}}). The layers are combined based on the \code{index}.

The function used should return a single value, and the number of layers in the output SpatRaster equals the number of unique values in \code{index}.

For example, if you have a SpatRaster with 6 layers, you can use \code{index=c(1,1,1,2,2,2)} and \code{fun=sum}. This will return a SpatRaster with two layers. The first layer is the sum of the first three layers in the input SpatRaster, and the second layer is the sum of the last three layers in the input SpatRaster. Indexes are recycled such that \code{index=c(1,2)} would also return a SpatRaster with two layers (one based on the odd layers (1,3,5), the other based on the even layers (2,4,6)).

See \code{\link{app}} or \code{\link{Summary-methods}} if you want to use a more efficient function that returns multiple layers based on \bold{all} layers in the SpatRaster object.
}

\usage{
\S4method{tapp}{SpatRaster}(x, index, fun, ..., filename="", overwrite=FALSE, wopt=list())
}

\arguments{
  \item{x}{SpatRaster}
  \item{index}{factor or numeric (integer). Vector of length \code{nlyr(x)} (shorter vectors are recycled) grouping the input layers}
  \item{fun}{function to be applied}
  \item{...}{additional arguments passed to \code{fun}}
  \item{filename}{character. Output filename}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list with named options for writing files as in \code{\link{writeRaster}}}
}

\value{
SpatRaster
}

\seealso{\code{\link{app}}, \code{\link{Summary-methods}}}

\examples{
r <- rast(ncol=10, nrow=10)
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

