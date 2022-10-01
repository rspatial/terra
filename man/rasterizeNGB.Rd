\name{rasterizeNGB}

\docType{methods}

\alias{rasterizeNGB}
\alias{rasterizeNGB,SpatRaster,SpatVector-method}
\alias{rasterizeNGB,SpatRaster,matrix-method}

\title{Rasterize points with a neighborhood}

\description{
Rasterize points using a ellipse (or circle) as neighborhood. For each raster cell, the points that fall within a user-defined ellipse that is centered on the cell are considered. 

This can result in much smoother results compared to the standard \code{\link{rasterize}} method.
}

\usage{
\S4method{rasterizeNGB}{SpatRaster,SpatVector}(x, y, field, algo, rad1, rad2=rad1, angle=0, minPoints=1, fill=NA, filename="", ...)

\S4method{rasterizeNGB}{SpatRaster,matrix}(x, y, algo, rad1, rad2=rad1, angle=0, minPoints=1, fill=NA, filename="", ...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{y}{SpatVector or matrix with three columns (x,y,z)}
  \item{field}{character. field name in SpatVector \code{y}}
  \item{algo}{character. The algorithm to use. One of: "min", "max", "mean", "range", "count", "distto", "distbetween". See Details}
  \item{rad1}{numeric. The first radius of the search ellipse (x-axis if \code{angle=0}) }
  \item{rad2}{numeric. The second radius of the search ellipse (y-axis if \code{angle=0}) }
  \item{angle}{numeric. The angle of ellipse rotation in degrees. Rotation is counter clockwise}
  \item{minPoints}{numeric. The minimum number of points to use. If fewer points are found in a search ellipse it is considered empty and the fill value is returned}
  \item{fill}{numeric. value to use to fill cells with empty search areas}
  \item{filename}{character. Output filename}
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
}


\seealso{
\code{\link{rasterize}}
}

\details{

Description of the algorithms

\tabular{ll}{\bold{Name} \tab \bold{Description} \cr
min    \tab Minimum value of z inside the ellipse\cr
max    \tab Maximum value or z inside the ellipse\cr
mean    \tab Mean value or z inside the ellipse\cr
range  \tab Range (max - min) of the values of z inside the ellipse\cr
count  \tab Number of points inside the ellipse  \cr
distto    \tab Average distance to the points inside the ellipse \cr
distbetween \tab  Average distance between points inside the ellipse\cr
}

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

r <- rasterizeNGB(r, xyz, "count", rad1=5)

rfuns <- c("count", "min", "max", "mean")
x <- lapply(rfuns, function(f) rasterizeNGB(r, xyz, f, rad1=5))
names(x) <- rfuns 
x <- rast(x)
plot(x)
}

\keyword{spatial}

