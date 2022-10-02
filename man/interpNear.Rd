\name{interpNear}

\docType{methods}

\alias{interpNear}
\alias{interpNear,SpatRaster,SpatVector-method}
\alias{interpNear,SpatRaster,matrix-method}

\title{Interpolate points using neighborhoods around raster cells}

\description{

near:  Nearest Neighbor

linear: Linear interpolation
}

\usage{
\S4method{interpNear}{SpatRaster,SpatVector}(x, y, field, algo, radius, fill=NA, filename="", ...)

\S4method{interpNear}{SpatRaster,matrix}(x, y, algo, radius, fill=NA, filename="", ...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{y}{SpatVector or matrix with three columns (x,y,z)}
  \item{field}{character. field name in SpatVector \code{y}}
  \item{algo}{character. The algorithm to use. One of: "nearest", "linear"}
  \item{radius}{numeric. The radius of the circle (single number). For \code{algo=nearest} only: if two numbers are supplied these are interpreted as the radii of an ellipse (x and y-axis). If a third number is supplied it is assumed to be in degrees and used to rotate the ellipse (counter clockwise)}

  \item{pars}{numeric. algorithm specific parameters. See Details}
  \item{fill}{numeric. value to use to fill empty cells}
  \item{filename}{character. Output filename}
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
}

\seealso{
\code{\link{rasterizeNGB}, \link{rasterize}, \link{interpIDW}, \link{interpolate}}
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

r <- interpolateNGB(r, xyz, "count", pars=c(Radius1=5,Radius2=5,Angle=0,MinPoints=0))

}

\keyword{spatial}

