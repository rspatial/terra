\name{zonal}

\alias{zonal}

\alias{zonal,SpatRaster,SpatRaster-method}

\title{Zonal statistics}

\description{
Compute zonal statistics, that is summarized values of a SpatRaster for each "zone" defined by another SpatRaster. 

If \code{fun} is a true \code{function}, \code{zonal} may fail for very large SpatRaster objects, except for the functions ("mean", "min", "max", or "sum"). 
}

\usage{
\S4method{zonal}{SpatRaster,SpatRaster}(x, z, fun=mean, ..., as.raster=FALSE, filename="", wopt=list()) 
}

\arguments{
  \item{x}{SpatRaster}
  \item{z}{SpatRaster with values representing zones}
  \item{fun}{function to be applied to summarize the values by zone. Either as character: "mean", "min", "max", "sum", or, for relatively small SpatRasters, a proper function}
  \item{...}{additional arguments passed to fun}  
  \item{as.raster}{logical. If \code{TRUE}, a SpatRaster is returned with the zonal statistic for each zone}  
  \item{filename}{character. Output filename (ignored if \code{as.raster=FALSE}}
  \item{wopt}{list with additional arguments for writing files as in \code{\link{writeRaster}}}
}

\value{
A \code{data.frame} with a value for each zone (unique value in \code{zones})
}

\seealso{ See \code{\link{global}} for "global" statistics (i.e., all of \code{x} is considered a single zone),  \code{\link{app}} for local statistics, and \code{\link{extract}} for summarizing values for polygons}

\examples{
r <- rast(ncols=10, nrows=10)
values(r) <- 1:ncell(r)
z <- rast(r)
values(z) <- rep(c(1:2, NA, 3:4), each=20)
names(z) <- "zone"
zonal(r, z, "sum", na.rm=TRUE)

# multiple layers
r <- rast(system.file("ex/logo.tif", package = "terra")) 
# zonal layer 
z <- rast(r, 1)
names(z) <- "zone"
values(z) <- rep(c(1:2, NA, c(3:4)), each=ncell(r)/5, length.out=ncell(r))

zonal(r, z, "mean", na.rm = TRUE)

# raster of zonal values
zr <- zonal(r, z, "mean", na.rm = TRUE, as.raster=TRUE)
}

\keyword{spatial}

