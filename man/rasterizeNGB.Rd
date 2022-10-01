\name{rasterizeNGB}

\docType{methods}

\alias{rasterizeNGB}
\alias{rasterizeNGB,SpatRaster,SpatVector-method}
\alias{rasterizeNGB,SpatRaster,matrix-method}

\title{Rasterize points using neighborhoods around raster cells}

\description{
}

\usage{
\S4method{rasterizeNGB}{SpatRaster,SpatVector}(x, y, field, algo, pars, filename="", ...)

\S4method{rasterizeNGB}{SpatRaster,matrix}(x, y, algo, pars, filename="", ...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{y}{SpatVector or matrix with three columns (x,y,z)}
  \item{field}{character. field name in SpatVector \code{y}}
  \item{algo}{character. The algorithm to use. One of: "min", "max", "range", "count", "distance", "ptsdistance", "invdistpow", "invdistpownear, "movingavg", "near", "linear". See Details}
  \item{pars}{numeric. algorithm specific parameters. See Details}
  \item{filename}{character. Output filename}
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
}

\seealso{
\code{\link{rasterize}}
}

\details{

explain algorithms and parameters 

"min", "max", "range", "count", "distance", "ptsdistance", "invdistpow", "invdistpownear, "movingavg", "near", "linear"

}

\value{
SpatRaster
}

\examples{
f <- system.file("ex/lux.shp", package="terra")
}

\keyword{spatial}

