\name{quantile}

\docType{methods}

\alias{quantile}
\alias{quantile,SpatRaster-method}
\alias{quantile,SpatVector-method}

\title{Quantiles of spatial data}

\description{
Compute quantiles for each cell across the layers of a SpatRaster. 

You can use use \code{\link{global}(x, fun=quantile)} to instead compute quantiles across cells for each layer.

You can also use this method to compute quantiles of the numeric variables of a SpatVector.
}

\usage{
\S4method{quantile}{SpatRaster}(x, probs=seq(0, 1, 0.25), na.rm=FALSE, filename="", ...) 

\S4method{quantile}{SpatVector}(x, probs=seq(0, 1, 0.25), ...) 
}

\arguments{
  \item{x}{SpatRaster or SpatVector}
  \item{probs}{numeric vector of probabilities with values in [0,1]}  
  \item{na.rm}{logical. If \code{TRUE}, \code{NA}'s are removed from \code{x} before the quantiles are computed}
  \item{filename}{character. Output filename}
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}  
}

\value{
SpatRaster with layers representing quantiles
}

\seealso{ \code{\link{app}}  }


\examples{
r <- rast(system.file("ex/logo.tif", package="terra"))   
rr <- c(r/2, r, r*2)
qr <- quantile(rr)
qr

\dontrun{
# same but slower
qa <- app(rr, quantile)
}

#quantile by layer instead of by cell
qg <- global(r, quantile)
}

\keyword{methods}
\keyword{spatial}

