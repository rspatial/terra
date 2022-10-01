\name{rasterizeNGB2}

\docType{methods}

\alias{rasterizeNGB2}
\alias{rasterizeNGB2,SpatRaster,SpatVector-method}
\alias{rasterizeNGB2,SpatRaster,matrix-method}

\title{Rasterize points using neighborhoods around raster cells}

\description{
}

\usage{
\S4method{rasterizeNGB2}{SpatRaster,SpatVector}(x, y, field, algo, pars, fill=NA, filename="", ...)

\S4method{rasterizeNGB2}{SpatRaster,matrix}(x, y, algo, pars, fill=NA, filename="", ...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{y}{SpatVector or matrix with three columns (x,y,z)}
  \item{field}{character. field name in SpatVector \code{y}}
  \item{algo}{character. The algorithm to use. One of: "invdistpow", "invdistpownear, "near", "linear". See Details}
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
invdistpow  \tab Inverse distance to a power \cr
invdistpownear \tab  Inverse distance to a power with nearest neighbor search for max points \cr
near        \tab Nearest Neighbor  \cr
linear      \tab  Linear interpolation (from Delaunay triangulation) \cr
}

Parameters required by algorithm

\tabular{lrrrrrr}{
\bold{Parameter}  \tab \bold{invdistpow} \tab \bold{invdistpownear} \tab \bold{near} \tab \bold{linear} \cr

Power      \tab x \tab x  \tab   \tab   \cr
Smoothing  \tab x \tab x  \tab   \tab   \cr
Radius     \tab   \tab x  \tab   \tab x \cr
Radius1    \tab x \tab    \tab x \tab   \cr
Radius2    \tab x \tab    \tab x \tab   \cr
Angle      \tab x \tab    \tab x \tab   \cr
nMaxPoints \tab x \tab x  \tab   \tab   \cr
nMinPoints \tab x \tab x  \tab   \tab   \cr
}

Description of the parameters 

\tabular{ll}{\bold{Parameter} \tab \bold{Description} \cr
Power     \tab Weighting power\cr
Smoothing \tab Smoothing parameter\cr
Radius    \tab The radius of the search circle \cr
Radius1   \tab The first radius (X axis if rotation angle is 0) of the search ellipse.\cr
Radius2   \tab The second radius (Y axis if rotation angle is 0) of the search ellipse.\cr
Angle     \tab Angle of ellipse rotation in degrees. Rotation is counter clockwise.\cr
MaxPoints \tab The maximum number of points to consider in a search area. Additional points are ignored. If fewer points are found, the fill value is returned\cr
MinPoints \tab The minimum number of points to use. If fewer points are found in a search area, the fill value is returned.\cr
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

r <- rasterizeNGB2(r, xyz, "count", pars=c(Radius1=5,Radius2=5,Angle=0,MinPoints=0))

}

\keyword{spatial}

