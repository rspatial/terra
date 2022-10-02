\name{rasterizeWin}

\docType{methods}

\alias{rasterizeWin}
\alias{rasterizeWin,SpatRaster,SpatVector-method}
\alias{rasterizeWin,SpatRaster,matrix-method}

\title{Rasterize points with a moving window}

\description{
Rasterize points using a circle (or ellipse) as moving window. For each raster cell, the points (\code{x, y}) that fall within the window centered on that cell are considered. A function is used to compute a summary value (e.g. "mean") for the values (\code{z}) associated with these points.

This can result in much smoother results compared to the standard \code{\link{rasterize}} method.
}

\usage{
\S4method{rasterizeWin}{SpatRaster,SpatVector}(x, y, field, fun, radius, minPoints=1, fill=NA, filename="", ...)

\S4method{rasterizeWin}{SpatRaster,matrix}(x, y, fun, radius, minPoints=1, fill=NA, filename="", ...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{y}{SpatVector or matrix with three columns (x, y, z)}
  \item{field}{character. field name in SpatVector \code{y} with the values to rasterize}
  \item{fun}{function to summarize the values for each cell. Functions must take a numeric vector and return a single number. The fast built-in functions are "min", "max", "mean", "range", "count", "distto" (average distance to the points), "distbetween" (average distance between the points)}
  \item{radius}{numeric. The radius of the circle (single number). If two numbers are supplied these are interpreted as the radii of an ellipse (x- and y-axis). If a third number is supplied it is assumed to be in degrees and used to rotate the ellipse (counter clockwise)}
  \item{minPoints}{numeric. The minimum number of points to use. If fewer points are found in a search ellipse it is considered empty and the fill value is returned}
  \item{fill}{numeric. value to use to fill cells with empty search areas}
  \item{filename}{character. Output filename}
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
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
xyz <- cbind(x,y,z)

r <- rasterizeWin(r, xyz, "count", radius=5)

rfuns <- c("count", "min", "max", "mean")
x <- lapply(rfuns, function(f) rasterizeWin(r, xyz, f, radius=5))
names(x) <- rfuns 
x <- rast(x)
#plot(x)


}

\keyword{spatial}

