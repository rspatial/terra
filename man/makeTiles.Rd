\name{makeTiles}

\docType{methods}


\alias{getTileExtents}
\alias{getTileExtents,SpatRaster-method}

\alias{makeTiles}
\alias{makeTiles,SpatRaster-method}


\title{Make tiles or get their extents}

\description{ 
Divide a SpatRaster into "tiles". The cells of another SpatRaster (normally with a much lower resolution) or a SpatVector with polygon geometry can be used to define the tiles. You can also provide one or two numbers to indicate the number of rows and columns per tile. 

\code{getTileExtents} returns the extents of the (virtual) tiles, while \code{makeTiles} creates files for the tiles and returns their filenames. 
}

\usage{
\S4method{makeTiles}{SpatRaster}(x, y, filename="tile_.tif", extend=FALSE,
		na.rm=FALSE, buffer=0, value="files", overwrite=FALSE, ...)


\S4method{getTileExtents}{SpatRaster}(x, y, extend=FALSE, buffer=0, cores=1)
}

\arguments{
  \item{x}{SpatRaster}
  \item{y}{SpatRaster or SpatVector defining the zones; or a positive integer specifying the number of rows and columns for each zone (or 2 numbers to differentiate the number of rows and columns). For \code{getTileExtents}, \code{y} may also be missing or \code{NULL}, in which case a tile size is computed automatically (see \emph{Details})}
  \item{filename}{character. Output filename template. Filenames will be altered by adding the tile number for each tile}
  \item{extend}{logical. If \code{TRUE}, the extent of \code{y} is expanded to assure that it covers all of \code{x}}
  \item{na.rm}{logical. If \code{TRUE}, tiles with only missing values are ignored}
  \item{buffer}{integer. The number of additional rows and columns added to each tile. Can be a single number, or two numbers to specify a separate number of rows and columns. This allows for creating overlapping tiles that can be used for computing spatial context dependent values with e.g. \code{\link{focal}}. The expansion is only inside \code{x}, no rows or columns outside of \code{x} are added}
  \item{cores}{integer (only used by \code{getTileExtents} when \code{y} is missing). The number of worker processes that will consume the tiles. Larger values lead to smaller per-worker tiles so that the concurrent peak memory stays within budget. See \emph{Details}}
  \item{value}{character. The type of return value desired. Either "files" (for the filenames), "raster" (for a SpatRaster), or "collection" (for a SpatRasterCollection)}
  \item{overwrite}{logical. If \code{TRUE}, existing tiles are overwritten; otherwise they are skipped (without error or warning)}
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
}

\details{
When \code{y} is missing in \code{getTileExtents}, the tile size is chosen so that

\itemize{
  \item tiles align to whole GDAL blocks of the source file(s) (as reported by \code{\link{fileBlocksize}}). When no source block size is reported, a default of 256 rows by 256 columns is used.
  \item the per-worker peak memory (\code{cells per tile} \eqn{\times} \code{nlyr} \eqn{\times} 8 bytes \eqn{\times} a small \code{ncopies} factor) stays within \code{cores} workers' share of \code{terraOptions("memfrac")} of free RAM (see \code{\link{free_RAM}}, \code{\link{terraOptions}}).
  \item there are at least a few tiles per worker, so the work load-balances.
  \item the result never exceeds the raster's own dimensions.
}

In practice this means that for a tiled GeoTIFF / COG, tiles are an integer multiple of the file block size; for an in-memory or non-tiled raster the result is roughly square.
}

\value{
\code{makeTiles} returns a character (filenames), SpatRaster or SpatRasterCollection value. \code{getTileExtents} returns a matrix with extents
}


\seealso{
\code{\link{vrt}} to create a SpatRaster from tiles; \code{\link{crop}} for sub-setting arbitrary parts of a SpatRaster; \code{\link{divide}} to divide a SpatRaster into parts.
}



\examples{
r <- rast(ncols=100, nrows=100)
values(r) <- 1:ncell(r)
x <- rast(ncols=2, nrows=2)

getTileExtents(r, x)
getTileExtents(r, x, buffer=3)

# auto: tile size from GDAL block size and a per-worker memory budget
getTileExtents(r)            # cores = 1
getTileExtents(r, cores=4)   # smaller tiles, sized for 4 concurrent workers


filename <- paste0(tempfile(), "_.tif")
ff <- makeTiles(r, x, filename)
ff

vrt(ff)
}


\keyword{methods}
\keyword{spatial}
