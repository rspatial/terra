\name{zonal}

\alias{zonal}

\alias{zonal,SpatRaster,SpatRaster-method}

\alias{zonal,SpatVector,SpatVector-method}

\title{Zonal statistics}

\description{
Compute zonal statistics, that is summarize values of a SpatRaster for each "zone" defined by another SpatRaster. 

If \code{fun} is a true R \code{function}, \code{zonal} may fail for very large SpatRaster objects, except for the functions ("mean", "min", "max", "sum", "isNA", and "notNA"). 

You can also summarize values of a SpatVector for each polygon (zone) defined by another SpatVector. 

See \code{\link{extract}} to summarize values of a SpatRaster with a SpatVector. 
}

\usage{
\S4method{zonal}{SpatRaster,SpatRaster}(x, z, fun=mean, ..., w=NULL, as.raster=FALSE, filename="", wopt=list()) 

\S4method{zonal}{SpatVector,SpatVector}(x, z, fun=mean, ..., weighted=FALSE, as.polygons=FALSE) 
}

\arguments{
  \item{x}{SpatRaster or SpatVector}
  \item{z}{Object of the same class as \code{x}. That is, either a SpatRaster with cell-values representing zones or a SpatVector with each polygon geometry representing a zone}
  \item{fun}{function to be applied to summarize the values by zone. Either as character: "mean", "min", "max", "sum", "isNA", and "notNA" and, for relatively small SpatRasters, a proper function}
  \item{...}{additional arguments passed to fun}  
  \item{w}{SpatRaster with weights. Should have a single-layer with non-negative values}
  \item{as.raster}{logical. If \code{TRUE}, a SpatRaster is returned with the zonal statistic for each zone}  
  \item{filename}{character. Output filename (ignored if \code{as.raster=FALSE}}
  \item{wopt}{list with additional arguments for writing files as in \code{\link{writeRaster}}}

  \item{weighted}{logical. If \code{TRUE}, a weighted.mean is computed and \code{fun} is ignored. Weights are based on the length of the lines or the area of the polygons in \code{x} that intersect with \code{z}. This argument is ignored of \code{x} is a SpatVector or points}  
  \item{as.polygons}{logical. Should the results be merged with the attributes of \code{z}?}
}

\value{
A \code{data.frame} with a value for each zone, or a SpatRaster, or SpatVector of polygons
}

\seealso{ See \code{\link{global}} for "global" statistics (i.e., all of \code{x} is considered a single zone),  \code{\link{app}} for local statistics, and \code{\link{extract}} for summarizing values for polygons}

\examples{
r <- rast(ncols=10, nrows=10)
values(r) <- 1:ncell(r)
z <- rast(r)
values(z) <- rep(c(1:2, NA, 3:4), each=20)
names(z) <- "zone"
zonal(r, z, "sum", na.rm=TRUE)

# with weights 
w <- init(r, "col")
zonal(r, z, w=w, "mean", na.rm=TRUE)


# multiple layers
r <- rast(system.file("ex/logo.tif", package = "terra")) 
# zonal layer 
z <- rast(r, 1)
names(z) <- "zone"
values(z) <- rep(c(1:2, NA, c(3:4)), each=ncell(r)/5, length.out=ncell(r))

zonal(r, z, "mean", na.rm = TRUE)

# raster of zonal values
zr <- zonal(r, z, "mean", na.rm = TRUE, as.raster=TRUE)


### SpatVector 

f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)[,c(2,4)]

p <- spatSample(v, 100)
values(p) <- data.frame(b2=1:100, ssep1=100:1)

zonal(p, v, mean)
}

\keyword{spatial}

