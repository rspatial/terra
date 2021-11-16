\name{scale}

\alias{scale}
\alias{scale,SpatRaster-method}


\title{Scale values}

\description{
Center and/or scale raster data. For details see \code{\link{scale}} 
}

\usage{
\S4method{scale}{SpatRaster}(x, center=TRUE, scale=TRUE)
}


\arguments{
 \item{x}{SpatRaster}
 \item{center}{logical or numeric. If \code{TRUE}, centering is done by subtracting the layer means (omitting \code{NA}s), and if \code{FALSE}, no centering is done. If \code{center} is a numeric vector (recycled to \code{nlyr(x)}), then each layer of \code{x} has the corresponding value from center subtracted from it.}
  \item{scale}{logical or numeric. If \code{TRUE}, scaling is done by dividing the (centered) layers of \code{x} by their standard deviations if \code{center} is \code{TRUE}, and the root mean square otherwise. If scale is \code{FALSE}, no scaling is done. If \code{scale} is a numeric vector (recycled to \code{nlyr(x)}), each layer of \code{x} is divided by the corresponding value. Scaling is done after centering.}
}
 
\value{
SpatRaster
}


\seealso{ \code{\link[base]{scale}} }

\examples{
r <- rast(system.file("ex/logo.tif", package="terra"))   
s <- scale(r)

## the equivalent, computed in steps
m <- global(r, "mean")
rr <- r - m[,1]
rms <- global(rr, "rms")
ss <- rr / rms[,1]

}

\keyword{ spatial }
