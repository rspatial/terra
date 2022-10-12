\name{sapp}

\docType{methods}

\alias{sapp}
\alias{sapp,SpatRaster-method}
\alias{sapp,SpatRasterDataset-method}

\title{Apply a terra function that takes only a single layer and returns a SpatRaster to all layers of a SpatRaster}

\description{
Apply to all layers of a SpatRaster a function that only takes a single layer SpatRaster and returns a SpatRaster (these are rare). In most cases you can also use \code{lapply} or \code{sapply} for this. 

Or apply the same method to each sub-dataset (SpatRaster) in a SpatRasterDataset
}

\usage{
\S4method{sapp}{SpatRaster}(x, fun, ..., filename="", overwrite=FALSE, wopt=list())

\S4method{sapp}{SpatRasterDataset}(x, fun, ..., filename="", overwrite=FALSE, wopt=list())
}

\arguments{
  \item{x}{SpatRaster or SpatRasterDataset}
  \item{fun}{if \code{x} is a \code{SpatRaster}: a function that takes a SpatRaster argument and can be applied to each layer of \code{x} (e.g. \code{\link{terrain}}. 
  if \code{x} is a \code{SpatRasterDataset}: a function that is applied to all layers of the SpatRasters in \code{x} (e.g. \code{mean}}
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

#SpatRasterDataset
sd <- sds(s*2, s/2)
y <- sapp(sd, mean)
z <- sapp(sd, \(i) 2 * mean(i))

}
 
\keyword{methods}
\keyword{spatial}

