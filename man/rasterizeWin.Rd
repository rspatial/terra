\name{rasterizeWin}

\docType{methods}

\alias{rasterizeWin}
\alias{rasterizeWin,SpatVector,SpatRaster-method}
\alias{rasterizeWin,data.frame,SpatRaster-method}

\title{Rasterize points with a moving window}

\description{
Rasterize points using a circle (or ellipse) as moving window. For each raster cell, the points (\code{x, y}) that fall within the window centered on that cell are considered. A function is used to compute a summary value (e.g. "mean") for the values (\code{z}) associated with these points.

This can result in much smoother results compared to the standard \code{\link{rasterize}} method.
}

\usage{
\S4method{rasterizeWin}{SpatVector,SpatRaster}(x, y, field, win="circle", pars, fun, ..., cvars=FALSE, 
			  minPoints=1, fill=NA, filename="", wopt=list())

\S4method{rasterizeWin}{data.frame,SpatRaster}(x, y, win="circle", pars, fun, ..., cvars=FALSE, 
		      minPoints=1, fill=NA, filename="", wopt=list())
}

\arguments{
  \item{x}{SpatVector or matrix with at least three columns ((x, y) coordinates and a variable to be rasterized)}
  \item{y}{SpatRaster}
  \item{field}{character. field name in SpatVector \code{x} with the values to rasterize}
  \item{win}{character to choose the window type. Can be "circle", "ellipse", "rectangle", or "buffer"}
  \item{pars}{parameters to define the window. If \code{win="circle"} or \code{win="buffer"}, a single number to set the radius of the circle or the width of the buffer. If \code{win="ellipse"}, either two numbers (the x and y-axis) or three numbers the axes and a rotation (in degrees). If \code{win="rectangle"}, either two (width, height) or three (width, height) and the rotation in degrees. The unit of the radius/width/height/axis parameters is that of the coordinate reference system (it is not expressed as cells). That is, if you have a lon/lat crs, there is no conversion of degrees to meters or vice-versa.}
  \item{fun}{function to summarize the values for each cell. If \code{cvars=FALSE}, functions must take a numeric vector and return (in all cases) one or more numbers. If \code{cvars=TRUE}, and multiple variables are used, the function must take a single argument (a data.frame with the names variables).   
  For \code{win="circle"} and \code{win="ellipse"} there are two additional character values that can be used: \code{"distto"} (average distance to the points from the center of the cell) and \code{"distbetween"} (average distance between the points inside the window)}
  \item{...}{additional named arguments passed to \code{fun}}
  \item{minPoints}{numeric. The minimum number of points to use. If fewer points are found in a search ellipse it is considered empty and the fill value is returned}
  \item{fill}{numeric. value to use to fill cells with empty search areas}
  \item{cvars}{logical. When using multiple fields, should \code{fun} operate on all of them at once? If not, \code{fun} is applied to each variable separately}
  \item{filename}{character. Output filename}
  \item{wopt}{list with additional arguments for writing files as in \code{\link{writeRaster}}}
}


\seealso{
\code{\link{rasterize}}, \code{\link{rasterizeGeom}}, \code{\link{interpNear}}, \code{\link{interpIDW}}
}


\value{
SpatRaster
}

\examples{

r <- rast(ncol=100, nrow=100, crs="local", xmin=0, xmax=50, ymin=0, ymax=50)
set.seed(100)
x <- runif(50, 5, 45)
y <- runif(50, 5, 45)
z <- sample(50)
xyz <- data.frame(x,y,z)

r <- rasterizeWin(xyz, r, fun="count", pars=5)

rfuns <- c("count", "min", "max", "mean")
x <- lapply(rfuns, function(f) rasterizeWin(xyz, r, fun=f, pars=5))
names(x) <- rfuns 
x <- rast(x)
#plot(x)

}

\keyword{spatial}


