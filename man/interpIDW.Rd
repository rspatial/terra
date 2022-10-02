\name{interpIDW}

\docType{methods}

\alias{interpIDW}
\alias{interpIDW,SpatRaster,SpatVector-method}
\alias{interpIDW,SpatRaster,matrix-method}

\title{Interpolate points using neighborhoods around raster cells}

\description{

invdistpow: Inverse distance to a power

invdistpownear: Inverse distance to a power with nearest neighbor search for max points

}

\usage{
\S4method{interpIDW}{SpatRaster,SpatVector}(x, y, field, algo, power, smooth, radius, minPoints, maxPoints, fill=NA, filename="", ...)

\S4method{interpIDW}{SpatRaster,SpatVector}(x, y, algo, power, smooth, radius, minPoints, maxPoints, fill=NA, filename="", ...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{y}{SpatVector or matrix with three columns (x,y,z)}
  \item{field}{character. field name in SpatVector \code{y}}
  \item{algo}{character. The algorithm to use. One of: "invdistpow", "invdistpownear}
  \item{radius}{numeric. The radius of the circle (single number). For \code{algo=invdistpow} only: if two numbers are supplied these are interpreted as the radii of an ellipse (x and y-axis). If a third number is supplied it is assumed to be in degrees and used to rotate the ellipse (counter clockwise)}
  \item{power}{numeric. Weighting power}
  \item{smooth}{numeric. Smoothing parameter}
  \item{minPoints}{numeric. The minimum number of points to use. If fewer points are found in a search ellipse it is considered empty and the fill value is returned}
  \item{maxPoints}{numeric. The maximum number of points to consider in a search area. Additional points are ignored. If fewer points are found, the fill value is returned}
  \item{fill}{numeric. value to use to fill empty cells}
  \item{filename}{character. Output filename}
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
}

\seealso{
\code{\link{rasterizeNGB}, \link{rasterize}, \link{interpNear}, \link{interpolate}}
}


\value{
SpatRaster
}

\examples{

r <- rast(ncol=100, nrow=100, crs="local", xmin=0, xmax=50, ymin=0, ymax=50)
set.seed(100)
x <- runif(250, 5, 45)
y <- runif(250, 5, 45)
z <- sample(250)
xyz <- cbind(x,y,z)

#r <- interpIDW(r, xyz, )

}

\keyword{spatial}

