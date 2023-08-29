\name{interpIDW}

\docType{methods}

\alias{interpIDW}
\alias{interpIDW,SpatRaster,SpatVector-method}
\alias{interpIDW,SpatRaster,matrix-method}

\title{Interpolate points using a moving window}

\description{

Interpolate points within a moving window using inverse distance weighting. The maximum number of points used can be restricted, optionally by selecting the nearest points. 
}

\usage{
\S4method{interpIDW}{SpatRaster,SpatVector}(x, y, field, radius, power=2, smooth=0,
       maxPoints=Inf, minPoints=1, near=TRUE, fill=NA, filename="", ...)

\S4method{interpIDW}{SpatRaster,matrix}(x, y, radius, power=2, smooth=0, 
       maxPoints=Inf, minPoints=1, near=TRUE, fill=NA, filename="", ...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{y}{SpatVector or matrix with three columns (x,y,z)}
  \item{field}{character. field name in SpatVector \code{y}}
  \item{radius}{numeric. The radius of the circle (single number). If \code{near=FALSE}, it is also possible to use two or three numbers. Two numbers are interpreted as the radii of an ellipse (x and y-axis). A third number should indicated the desired, counter clockwise, rotation of the ellipse (in degrees)}

  \item{power}{numeric. Weighting power}
  \item{smooth}{numeric. Smoothing parameter}
  \item{minPoints}{numeric. The minimum number of points to use. If fewer points are found in a search ellipse it is considered empty and the fill value is returned}
  \item{maxPoints}{numeric. The maximum number of points to consider in a search area. Additional points are ignored. If fewer points are found, the fill value is returned}
  \item{near}{logical. Should the nearest points within the neighborhood be used if \code{maxPoints} is reached?}
  \item{fill}{numeric. value to use to fill empty cells}
  \item{filename}{character. Output filename}
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
}

\seealso{
\code{\link{rasterizeWin}, \link{rasterize}, \link{interpNear}, \link{interpolate}}
}


\value{
SpatRaster
}

\examples{

r <- rast(ncol=100, nrow=100, crs="local", xmin=0, xmax=50, ymin=0, ymax=50)
set.seed(100)
x <- runif(25, 5, 45)
y <- runif(25, 5, 45)
z <- sample(25)
xyz <- cbind(x,y,z)

x <- interpIDW(r, xyz, radius=5, power=1, smooth=1, maxPoints=5)

}

\keyword{spatial}

