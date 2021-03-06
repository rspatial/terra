\name{shade}

\alias{shade}

\title{Hill shading}

\description{
Compute hill shade from slope and aspect layers (both in radians). Slope and aspect can be computed with function \code{\link{terrain}}. 

A hill shade layer is often used as a backdrop on top of which another, semi-transparent, layer is drawn.
}

\usage{
shade(slope, aspect, angle=45, direction=0, normalize=FALSE, filename="", ...)  
}

\arguments{
  \item{slope}{SpatRasterwith slope values (in radians) }
  \item{aspect}{SpatRaster with aspect values (in radians) }
  \item{angle}{ The the elevation angle of the light source (sun), in degrees}
  \item{direction}{ The direction (azimuth) angle of the light source (sun), in degrees}
  \item{normalize}{Logical. If \code{TRUE}, values below zero are set to zero and the results are multiplied with 255}
  \item{filename}{character. Output filename}
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}  
}

\seealso{ \code{\link{terrain}} }

\references{
Horn, B.K.P., 1981. Hill shading and the reflectance map. Proceedings of the IEEE 69(1):14-47
}


\examples{
f <- system.file("ex/elev.tif", package="terra")
r <- rast(f)
alt <- disaggregate(r, 10, method="bilinear")
slope <- terrain(alt, 'slope')
aspect <- terrain(alt, 'aspect')
hill <- shade(slope, aspect, 40, 270)
plot(hill, col=grey(0:100/100), legend=FALSE, mar=c(2,2,1,4))
plot(alt, col=rainbow(25, alpha=0.35), add=TRUE)
}


\keyword{spatial}

