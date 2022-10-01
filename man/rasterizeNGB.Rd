\name{rasterizeNGB}

\docType{methods}

\alias{rasterizeNGB}
\alias{rasterizeNGB,SpatRaster,SpatVector-method}
\alias{rasterizeNGB,SpatRaster,matrix-method}

\title{Rasterize points using neighborhoods around raster cells}

\description{
}

\usage{
\S4method{rasterizeNGB}{SpatRaster,SpatVector}(x, y, field, algo, pars, fill=NA, filename="", ...)

\S4method{rasterizeNGB}{SpatRaster,matrix}(x, y, algo, pars, fill=NA, filename="", ...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{y}{SpatVector or matrix with three columns (x,y,z)}
  \item{field}{character. field name in SpatVector \code{y}}
  \item{algo}{character. The algorithm to use. One of: "min", "max", "range", "count", "distance", "ptsdistance", "invdistpow", "invdistpownear, "movingavg", "near", "linear". See Details}
  \item{pars}{numeric. algorithm specific parameters. See Details}
  \item{fill}{numeric. value to use to fill empty cells}
  \item{filename}{character. Output filename}
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
}

\seealso{
\code{\link{rasterize}}
}

\details{

Description of the algorithms

\tabular{ll}{\bold{Name} \tab \bold{Description} \cr
min    \tab Minimum Value  \cr
max    \tab Maximum Value  \cr
range  \tab Data Range  \cr
count  \tab Number of Points  \cr
distance    \tab Average Distance  \cr
ptsdistance \tab  Average Distance Between Data Points \cr
invdistpow  \tab Inverse distance to a power  \cr
invdistpownear \tab  Inverse distance to a power with nearest neighbor search for max points \cr
movingavg   \tab Moving Average  \cr
near        \tab Nearest Neighbor  \cr
linear      \tab  Linear interpolation (from Delaunay triangulation) \cr
}

Parameters required by algorithm

\tabular{lrrrrrr}{
\bold{Parameter} \tab \bold{min, max, range, count} \tab \bold{invdistpow} \tab \bold{invdistpownear} \tab \bold{movingavg} \tab \bold{near} \tab \bold{linear} \cr
          \tab \bold{distance, ptsdistance} \tab  \tab  \tab   \tab   \tab   \cr

Power      \tab   \tab x \tab x \tab   \tab   \tab   \cr
Smoothing  \tab   \tab x \tab x \tab   \tab   \tab   \cr
Radius     \tab   \tab   \tab x \tab   \tab   \tab x \cr
Radius1    \tab x \tab x \tab   \tab x \tab x \tab   \cr
Radius2    \tab x \tab x \tab   \tab x \tab x \tab   \cr
Angle      \tab x \tab x \tab   \tab x \tab x \tab x \cr
nMaxPoints \tab   \tab x \tab x \tab   \tab   \tab   \cr
nMinPoints \tab x \tab x \tab x \tab x \tab   \tab   \cr
}

Description of the parameters 

\tabular{ll}{\bold{Parameter} \tab \bold{Description} \cr
Power     \tab Weighting power\cr
Smoothing \tab Smoothing parameter\cr
Radius    \tab The radius of search circle \cr
Radius1   \tab The first radius (X axis if rotation angle is 0) of search ellipse.\cr
Radius2   \tab The second radius (Y axis if rotation angle is 0) of search ellipse.\cr
Angle     \tab Angle of ellipse rotation in degrees. Ellipse rotated counter clockwise.\cr
MaxPoints \tab Maximum number of data points to use. Do not search for more points than this number. If less amount of points found the grid node considered empty and will be filled with NODATA marker.\cr
MinPoints \tab Minimum number of data points to use. If less amount of points found the grid node considered empty and will be filled with NODATA marker.\cr
}

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

r <- rasterizeNGB(r, xyz, "count", pars=c(Radius1=5,Radius2=5,Angle=0,MinPoints=0))

}

\keyword{spatial}

