\name{shade}

\alias{shade}

\title{Hill shading}

\description{
Compute hill shade from slope and aspect layers (both in radians). Slope and aspect can be computed with function \code{\link{terrain}}. 

A hill shade layer is often used as a backdrop on top of which another, semi-transparent, layer is drawn.
}

\usage{
shade(slope, aspect, angle=45, direction=0, normalize=FALSE,
   filename="", overwrite=FALSE, ...)  
}

\arguments{
  \item{slope}{SpatRasterwith slope values (in radians) }
  \item{aspect}{SpatRaster with aspect values (in radians) }
  \item{angle}{ The elevation angle(s) of the light source (sun), in degrees}
  \item{direction}{ The direction (azimuth) angle(s) of the light source (sun), in degrees}
  \item{normalize}{Logical. If \code{TRUE}, values below zero are set to zero and the results are multiplied with 255}
  \item{filename}{character. Output filename}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}  
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}  
}

\seealso{ \code{\link{terrain}} }

\references{
Horn, B.K.P., 1981. Hill shading and the reflectance map. Proceedings of the IEEE 69(1):14-47
}


\examples{
f <- system.file("ex/elev.tif", package="terra")
r <- rast(f)
alt <- disagg(r, 10, method="bilinear")
slope <- terrain(alt, "slope", unit="radians")
aspect <- terrain(alt, "aspect", unit="radians")
hill <- shade(slope, aspect, 40, 270)
plot(hill, col=grey(0:100/100), legend=FALSE, mar=c(2,2,1,4))
plot(alt, col=rainbow(25, alpha=0.35), add=TRUE)

# A better hill shade may be achieved by combining 
# different angles and directions. For example

h <- shade(slope, aspect, angle = c(45, 45, 45, 80), direction = c(225, 270, 315, 135))
h <- Reduce(mean, h)				
}


\keyword{spatial}

