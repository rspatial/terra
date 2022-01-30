\name{mosaic}

\docType{methods}

\alias{mosaic}
\alias{mosaic,SpatRaster,SpatRaster-method}
\alias{mosaic,SpatRasterCollection,missing-method}

\title{ 
mosaic SpatRasters
}

\description{
Combine adjacent and (partly) overlapping SpatRasters to form a single new SpatRaster. Values in overlapping cells are averaged. 

This method is similar to the simpler, but faster \code{\link{merge}} method. 
}

\usage{
\S4method{mosaic}{SpatRaster,SpatRaster}(x, y, ..., fun="mean", filename="", overwrite=FALSE, wopt=list())

\S4method{mosaic}{SpatRasterCollection,missing}(x, fun="mean", filename="", ...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{y}{object of same class as \code{x}}
  \item{...}{additional SpatRasters}
  \item{fun}{character. One of "sum", "mean", "median", "min", "max"}
  \item{filename}{character. Output filename}
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
x <- rast(xmin=-110, xmax=-80, ymin=40, ymax=70, ncols=30, nrows=30)
y <- rast(xmin=-85, xmax=-55, ymax=60, ymin=30, ncols=30, nrows=30)
z <- rast(xmin=-60, xmax=-30, ymax=50, ymin=20, ncols=30, nrows=30)
values(x) <- 1:ncell(x)
values(y) <- 1:ncell(y)
values(z) <- 1:ncell(z)

m1 <- mosaic(x, y, z)
m2 <- mosaic(z, y, x)

# if you have many SpatRasters make a SpatRasterCollection from a list
rlist <- list(x, y, z)
rsrc <- sprc(rlist)

m <- mosaic(rsrc)
}

\keyword{methods}
\keyword{spatial}

