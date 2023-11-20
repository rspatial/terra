\name{zonal}

\alias{zonal}

\alias{zonal,SpatRaster,SpatRaster-method}
\alias{zonal,SpatRaster,SpatVector-method}
\alias{zonal,SpatVector,SpatVector-method}

\title{Zonal statistics}

\description{
Compute zonal statistics, that is summarize values of a SpatRaster for each "zone" defined by another SpatRaster, or by a SpatVector with polygon geometry. 

If \code{fun} is a true R \code{function}, the <SpatRaster,SpatRaster> method may fail when using very large SpatRasters, except for the functions ("mean", "min", "max", "sum", "isNA", and "notNA"). 

You can also summarize values of a SpatVector for each polygon (zone) defined by another SpatVector. 
}

\usage{
\S4method{zonal}{SpatRaster,SpatRaster}(x, z, fun="mean", ..., w=NULL, wide=TRUE,
		as.raster=FALSE, filename="", overwrite=FALSE, wopt=list())

\S4method{zonal}{SpatRaster,SpatVector}(x, z, fun="mean", na.rm=FALSE, w=NULL, weights=FALSE,
		exact=FALSE, touches=FALSE, as.raster=FALSE, as.polygons=FALSE, wide=TRUE, 
		filename="", wopt=list())

\S4method{zonal}{SpatVector,SpatVector}(x, z, fun=mean, ..., weighted=FALSE, as.polygons=FALSE) 
}

\arguments{
  \item{x}{SpatRaster or SpatVector}
  \item{z}{SpatRaster with cell-values representing zones or a SpatVector with each polygon geometry representing a zone. \code{z} can have multiple layers to define intersecting zones}
  \item{fun}{function to be applied to summarize the values by zone. Either as character: "mean", "min", "max", "sum", "isNA", and "notNA" and, for relatively small SpatRasters, a proper function}
  \item{...}{additional arguments passed to fun, such as \code{na.rm=TRUE}}  
  \item{w}{SpatRaster with weights. Should have a single-layer with non-negative values}
  \item{wide}{logical. Should the values returned in a wide format? For the \code{SpatRaster, SpatRaster} method this only affects the results when \code{nlyr(z) == 2}. For the \code{SpatRaster, SpatVector} method this only affects the results when \code{fun=table}}
  \item{as.raster}{logical. If \code{TRUE}, a SpatRaster is returned with the zonal statistic for each zone}  
  \item{filename}{character. Output filename (ignored if \code{as.raster=FALSE}}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}
  \item{wopt}{list with additional arguments for writing files as in \code{\link{writeRaster}}}

  \item{weights}{logical. If \code{TRUE} and \code{y} has polygons, the approximate fraction of each cell that is covered is returned as well, for example to compute a weighted mean}
  \item{exact}{logical. If \code{TRUE} and \code{y} has polygons, the exact fraction of each cell that is covered is returned as well, for example to compute a weighted mean}
  \item{touches}{logical. If \code{TRUE}, values for all cells touched by lines or polygons are extracted, not just those on the line render path, or whose center point is within the polygon. Not relevant for points; and always considered \code{TRUE} when \code{weights=TRUE} or \code{exact=TRUE}}

  \item{weighted}{logical. If \code{TRUE}, a weighted.mean is computed and \code{fun} is ignored. Weights are based on the length of the lines or the area of the polygons in \code{x} that intersect with \code{z}. This argument is ignored of \code{x} is a SpatVector or points}  
  \item{as.polygons}{logical. Should the zonal statistics be combined with the geometry of \code{z}?}
  \item{na.rm}{logical. If \code{TRUE}, \code{NA}s are removed}
}

\value{
A \code{data.frame} with a value for each zone, or a SpatRaster, or SpatVector of polygons.
}

\seealso{ See \code{\link{global}} for "global" statistics (i.e., all of \code{x} is considered a single zone), \code{\link{app}} for local statistics, and \code{\link{extract}} for an alternative way to summarize values of a SpatRaster with a SpatVector. With \code{\link{aggregate}} you can compute statistics for cell blocks defined by a number of rows and columns. 
}

\examples{

### SpatRaster, SpatRaster
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


### SpatRaster, SpatVector
x <- rast(ncol=2,nrow=2, vals=1:4, xmin=0, xmax=1, ymin=0, ymax=1, crs="+proj=utm +zone=1")
p <- as.polygons(x)
pp <- shift(p, .2)
r <- disagg(x, 4)

zonal(r, p)
zonal(r, p, sum)
zonal(x, pp, exact=TRUE)
zonal(c(x, x*10), pp, w=x)


### SpatVector, SpatVector

f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)[,c(2,4)]

p <- spatSample(v, 100)
values(p) <- data.frame(b2=1:100, ssep1=100:1)

zonal(p, v, mean)
}

\keyword{spatial}

