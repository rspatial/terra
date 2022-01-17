\name{sapp}

\docType{methods}

\alias{sapp}
\alias{sapp,SpatRaster-method}
\alias{sapp,SpatRasterDataset-method}

\title{Apply a terra function that takes only a single layer and returns a SpatRaster to all layers of a SpatRaster}

\description{
Apply to all layers of a SpatRaster a function that only takes a single layer SpatRaster and returns a SpatRaster (these are rare). 

In most cases you can also use \code{lapply} or \code{sapply} for this. 
}

\usage{
\S4method{sapp}{SpatRaster}(x, fun, ..., filename="", overwrite=FALSE, wopt=list())
}

\arguments{
  \item{x}{SpatRaster}
  \item{fun}{a function that takes a SpatRaster argument and can be applied to each layer of \code{x}}
  \item{...}{additional arguments to be passed to \code{fun}} 
  \item{filename}{character. Output filename}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list with named options for writing files as in \code{\link{writeRaster}}}
}

\value{
SpatRaster
}

\seealso{\code{ \link{lapp}, \link{app}, \link[terra]{tapp}, \link{lapply}} }


\examples{
s <- rast(system.file("ex/logo.tif", package="terra")) + 1  
x <- sapp(s, terrain)
}
 
\keyword{methods}
\keyword{spatial}

