\name{mosaic}

\docType{methods}

\alias{mosaic}
\alias{mosaic,SpatRaster,SpatRaster-method}

\title{ 
mosaic SpatRasters
}

\description{
mosaic SpatRasters to form a new SpatRaster object with a larger spatial extent. Unlike with \code{\link{merge}}, values on overlapping areas are averaged. 
}

\usage{
\S4method{mosaic}{SpatRaster,SpatRaster}(x, y, ..., fun="mean", filename="", overwrite=FALSE, wopt=list())
}

\arguments{
  \item{x}{SpatRaster}
  \item{y}{object of same class as \code{x}}
  \item{...}{additional SpatRasters}
  \item{fun}{character. One of "sum", "mean", "median", "min", "max"}
  \item{filename}{character. Output filename. Optional}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list with named options for writing files as in \code{\link{writeRaster}}}
} 

\details{
The SpatRaster objects must have the same origin and spatial resolution. In areas where the SpatRaster objects overlap, the values of the SpatRaster that is last in the sequence of arguments will be retained. 
}

\value{
SpatRaster
}

\seealso{
\code{\link{merge}}
}


\examples{
x <- rast(xmin=-110, xmax=-50, ymin=40, ymax=70, ncols=60, nrows=30)
y <- rast(xmin=-80, xmax=-20, ymax=60, ymin=30)
res(y) <- res(x)
values(x) <- 1:ncell(x)
values(y) <- 1:ncell(y)
mr <- mosaic(x, y)
}

\keyword{methods}
\keyword{spatial}

