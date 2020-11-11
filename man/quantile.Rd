\name{quantile}

\docType{methods}

\alias{quantile}
\alias{quantile,SpatRaster-method}

\title{SpatRaster local quantiles}

\description{
Compute quantiles for each cell across the layers of a SpatRaster 
}

\usage{
\S4method{quantile}{SpatRaster}(x, probs=seq(0, 1, 0.25), na.rm=FALSE,
    filename="", overwrite=FALSE, wopt=list(), ...) 
}

\arguments{
  \item{x}{SpatRaster}
  \item{probs}{numeric vector of probabilities with values in [0,1]}  
  \item{na.rm}{logical. If \code{TRUE}, \code{NA}'s are removed from \code{x} before the quantiles are computed}
  \item{filename}{character. Output filename. Optional}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list. Options for writing files as in \code{\link{writeRaster}}}  
  \item{...}{additional arguments, none implemented}
}

\value{
SpatRaster with layers representing quantiles
}

\seealso{ \code{\link{app}}  }



\examples{
r <- rast(system.file("ex/logo.tif", package="terra"))   
r <- c(r/2, r, r*2)
q <- quantile(r)
q

# same as (but faster than)
qa <- app(r, quantile)
}

\keyword{methods}
\keyword{spatial}

