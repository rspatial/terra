\name{update}

\alias{update}
\alias{update,SpatRaster-method}

  
\title{Update a raster file}

\description{
Update metadata or cell values of a file that is the data source of a SpatRaster. This modifies the file on disk directly without reading the entire raster into memory.

BE CAREFUL with this function as you are overwriting data in an existing file.
}

\usage{
\S4method{update}{SpatRaster}(object, crs=FALSE, extent=FALSE, names=FALSE,
    cells=NULL, values=NULL, layer=0)
}

\arguments{
  \item{object}{SpatRaster with a file source}
  \item{crs}{logical. Should the coordinate reference system be updated?}
  \item{extent}{logical. Should the extent be updated?}
  \item{names}{logical. Should the names be updated?}
  \item{cells}{numeric. Cell numbers to update (1-indexed). Must be used together with \code{values}}
  \item{values}{numeric. New values for the specified cells. Can be a single value (recycled), one value per cell (applied to all layers), or \code{length(cells) * length(layer)} values (one per cell per layer, in layer-major order)}
  \item{layer}{numeric or character. Layer(s) to update. Use \code{0} (default) to update all layers. Positive integers or layer names select specific layers}
}

\value{SpatRaster (invisibly)}

\details{
The \code{cells} and \code{values} arguments allow updating specific cell values on disk without loading the raster into memory. This is useful for modifying very large rasters that cannot fit in memory. See \code{\link{set.values}} for modifying values of in-memory rasters.
}

\examples{
## update metadata
s <- rast(system.file("ex/logo.tif", package="terra"))   
fname <- paste0(tempfile(), ".tif")
x <- writeRaster(s*1, fname)

ext(x) <- ext(x) + 1
crs(x) <- "+proj=utm +zone=1"
names(x) <- LETTERS[1:3]

update(x, crs=TRUE, extent=TRUE, names=TRUE)

rast(fname)

## update cell values on disk
r <- rast(nrows=10, ncols=10, vals=1:100)
f <- paste0(tempfile(), ".tif")
x <- writeRaster(r, f)
x[c(1, 2, 50, 51, 100)]
update(x, cells=c(1, 50, 100), values=c(999, 888, 777))
rast(f)[c(1, 2, 50, 51, 100)]

## update specific layer only
r <- rast(nrows=5, ncols=5, nlyr=3, vals=1:75)
f <- paste0(tempfile(), ".tif")
x <- writeRaster(r, f, datatype="FLT4S")
x[c(1, 2, 25, 26)]
update(x, cells=c(1, 25), values=c(0, 0), layer=2)
rast(f)[c(1, 2, 25, 26)]
}

\keyword{spatial}

