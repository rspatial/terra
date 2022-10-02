\name{interpNear}

\docType{methods}

\alias{interpNear}
\alias{interpNear,SpatRaster,SpatVector-method}
\alias{interpNear,SpatRaster,matrix-method}

\title{Nearest neighbor interpolation}

\description{
Nearest neighbor interpolation of points, using a moving window
}

\usage{
\S4method{interpNear}{SpatRaster,SpatVector}(x, y, field, radius, interpolate=FALSE, fill=NA, filename="", ...)

\S4method{interpNear}{SpatRaster,matrix}(x, y, radius, interpolate=FALSE, fill=NA, filename="", ...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{y}{SpatVector or matrix with three columns (x,y,z)}
  \item{field}{character. field name in SpatVector \code{y}}
  \item{radius}{numeric. The radius of the circle (single number). If \code{interpolate=FALSE} it is also possible to use two or three numbers. Two numbers are interpreted as the radii of an ellipse (x and y-axis). A third number should indicated the desired, counter clockwise, rotation of the ellipse (in degrees)}  
  \item{interpolate}{logical. Should the nearest neighbor values be linearly interpolated between points?}
  \item{fill}{numeric. value to use to fill empty cells}
  \item{filename}{character. Output filename}
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
}

\seealso{
\code{\link{rasterizeWin}, \link{rasterize}, \link{interpIDW}, \link{interpolate}}
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

x <- interpNear(r, xyz, radius=5)

p <- vect(data.frame(xyz), geom=c("x", "y")) 
v <- voronoi(p)

plot(x, col=rainbow(25))
lines(v)

# plot(v, col=rainbow(25)); points(p)
}

\keyword{spatial}

