\name{tile_apply}

\alias{tile_apply}

\title{Apply a function to tiles of a SpatRaster in parallel}

\description{
A high-level convenience wrapper that splits a \code{SpatRaster} into tiles, applies a user function to each tile (with \code{\link{window}} set under the hood so no values are copied), and assembles the per-tile results back into a single \code{SpatRaster}. The work can be distributed over a cluster of worker processes.

This is the \dQuote{Strategy 3} pattern from the \pkg{terra} parallelization guide, packaged as one call.
}

\usage{
tile_apply(x, fun, cores=1, cpkgs=NULL, tiles=NULL, buffer=0, ...,
           filename="", overwrite=FALSE, wopt=list(),
           overlap_fun=NULL)
}

\arguments{
  \item{x}{\code{SpatRaster}. Should normally be backed by a file (so that workers can re-open it cheaply); in-memory rasters are shipped to the workers via \code{\link{wrap}}}
  \item{fun}{function. Takes a \code{SpatRaster} as its first argument and must return a \code{SpatRaster}. It is called once per tile, with the windowed \code{x} as the first argument and any extra \code{...} arguments appended}
  \item{cores}{integer, or a \code{cluster} object created with \code{parallel::makeCluster}. If a number greater than 1, a PSOCK cluster of that size is created and torn down inside the call}
  \item{cpkgs}{character. Names of packages that must be loaded on each worker before \code{fun} is called (e.g. \code{"randomForest"})}
  \item{tiles}{specification of the tiles. If \code{NULL} (the default), the tile size is chosen automatically from the GDAL block size of \code{x} and a per-worker memory budget that depends on \code{cores} (see \code{\link{getTileExtents}} \emph{Details}). Otherwise one of: a \code{SpatExtent}; a list of \code{SpatExtent}s; a 4-column matrix of extents \code{(xmin, xmax, ymin, ymax)} as returned by \code{\link{getTileExtents}}; a \code{SpatRaster} or \code{SpatVector} that defines the tile geometry; or one or two integers specifying the number of rows and columns \emph{per tile} (passed on to \code{\link{getTileExtents}})}
  \item{buffer}{integer. Number of additional rows and columns (one number, or two for separate row/column buffers) of \code{x} to read on each side of every tile, so that \code{fun} can compute neighbourhood-dependent values (e.g. \code{\link{focal}}) without edge artefacts. The result of \code{fun} is then cropped back to the un-buffered tile extent before it is written, so the per-tile outputs are non-overlapping and assemble cleanly. Only used when \code{tiles = NULL}; ignored with a warning otherwise (when supplying \code{tiles} explicitly, build the overlap into the tiles yourself with \code{getTileExtents(..., buffer=)} and pick an \code{overlap_fun})}
  \item{...}{additional arguments passed to \code{fun}. When \code{cores > 1} these are serialized and shipped to each worker; \code{SpatRaster}, \code{SpatVector} and other \pkg{terra} objects in \code{...} are automatically wrapped and unwrapped}
  \item{filename}{character. Output filename for the assembled result. When empty (the default) and \code{overlap_fun = NULL}, a virtual raster (VRT) is returned that references the per-tile files (see \emph{Details})}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten if it exists}
  \item{wopt}{list. Writing options as in \code{\link{writeRaster}}, applied both to the per-tile intermediate files and to the final output}
  \item{overlap_fun}{character or \code{NULL}. The function name used by \code{\link{mosaic}} to combine overlapping cells (e.g. \code{"mean"}, \code{"min"}, \code{"max"}, \code{"sum"}, \code{"first"}). Use this when \code{tiles} was built with overlap (e.g. \code{getTileExtents(..., buffer=)}). When \code{NULL} (the default), tiles are assembled with \code{\link{vrt}}, which is correct for non-overlapping tiles (including the auto path with \code{buffer}, which crops away the overlap) and is essentially free}
}

\details{
For each tile, \code{tile_apply} sets a read window on \code{x} (\code{window(x) <- e}), so no pixel values are copied: subsequent operations on the windowed raster behave as if it had only that extent. The user function \code{fun} therefore sees a normal \code{SpatRaster} restricted to the tile.

\strong{Memory safety.} Each tile's result is written to its own file in a per-call sub-directory of \code{tempdir()} (as a GeoTIFF, honouring \code{wopt}). Workers therefore return only a filename, not pixel data, so the parent process never holds more than one tile's worth of values in RAM. The sequential path follows the same contract.

\strong{Buffered tiles for focal-style operations.} When \code{tiles = NULL} and \code{buffer > 0}, each tile is read on a slightly expanded extent (\code{buffer} cells on each side, clamped to \code{x}'s extent), \code{fun} is applied, and the result is cropped back to the un-buffered tile extent before it is written. This avoids the edge effects that operations like \code{\link{focal}} would otherwise produce at tile boundaries, while still letting the assembled output use the cheap, lossless \code{\link{vrt}} path. The buffer also enlarges the per-worker peak RAM (a tile reads \code{(rows + 2 buffer)(cols + 2 buffer)} cells), which is usually negligible since \code{buffer} is small compared with the auto tile dimensions, but if you choose a large buffer you may want to raise \code{cores} so the auto sizer picks smaller tiles.

The per-tile files are then assembled into the final \code{SpatRaster}:
\itemize{
  \item When \code{overlap_fun = NULL} (the default), the assembly is a virtual raster (\code{\link{vrt}}). This is essentially free and never reads the pixel data again. It is correct as long as the tiles do not overlap; if they do (e.g. when built with \code{getTileExtents(..., buffer=)}), GDAL keeps the value of the last tile drawn in the overlap region.
  \item When \code{overlap_fun} is set (e.g. \code{"mean"}), the assembly uses \code{\link{mosaic}} which streams the data through terra's standard chunked I/O and applies \code{overlap_fun} on overlapping cells. This is the right path for buffered tiles, focal operations, and similar.
}

If \code{filename} is supplied, the assembled result is materialised there and the per-tile intermediate files are removed. If \code{filename} is empty and a VRT is being built, the per-tile files are kept for the rest of the R session because they back the returned raster.

When \code{cores > 1}, \code{x} is shipped to the workers with \code{\link{wrap}}. For file-backed rasters this is cheap (only the metadata and filename travel across the worker boundary, and each worker re-opens the file independently). For in-memory rasters the values are copied to each worker.

Extra arguments passed via \code{...} must be named when \code{cores > 1}.
}

\value{
A \code{SpatRaster}.
}

\seealso{
\code{\link{getTileExtents}}, \code{\link{makeTiles}}, \code{\link{mosaic}}, \code{\link{vrt}}, \code{\link{window}}, \code{\link{wrap}}, \code{\link{app}}
}

\examples{
\donttest{
f <- system.file("ex/elev.tif", package="terra")
r <- rast(f)

# auto: block- and memory-aware tile size, sized for 1 worker.
# Returned object is a VRT backed by per-tile files in tempdir().
out1 <- tile_apply(r, function(x) x * 2)

# write the assembled result straight to disk (per-tile intermediates are
# removed automatically once the output file is materialised).
out2 <- tile_apply(r, function(x) x + 1,
                   filename=file.path(tempdir(), "out2.tif"), overwrite=TRUE)

# focal with the auto path: just ask for a buffer wide enough for the
# focal window. tile_apply reads (tile + buffer) cells per tile, runs the
# focal, crops back to the tile, and assembles via vrt() - no overlap_fun
# needed.
out3 <- tile_apply(r, function(x) focal(x, w=5, fun="mean", na.rm=TRUE),
                   buffer=2)

# focal with explicit pre-buffered tiles: use overlap_fun to blend overlap
tiles_buf <- getTileExtents(r, 50, buffer=3)
out3b <- tile_apply(r, function(x) focal(x, w=5, fun="mean", na.rm=TRUE),
                    tiles=tiles_buf, overlap_fun="mean")

# parallel on 2 workers (auto tiles are sized down for 2 concurrent workers)
\dontrun{
out4 <- tile_apply(r, function(x) x * 2, cores=2)

# pass extra arguments by name; nested terra objects are auto-wrapped
out5 <- tile_apply(r, function(x, k) x * k, cores=2, tiles=50, k=10)
}
}
}

\keyword{spatial}
