\name{tile_apply}

\alias{tile_apply}

\title{Parallelize a SpatRaster function}

\description{
The function splits a \code{SpatRaster} into virtual tiles and applies a user supplied function to each tile, and assembles the per-tile results back into a single \code{SpatRaster}. The work can be distributed over a cluster of worker processes.
}

\usage{
tile_apply(x, fun, cores=1, cpkgs=NULL, tiles=NULL, buffer=0, ...,
	filename="", overwrite=FALSE, wopt=list(), overlap_fun=NULL)
}

\arguments{
  \item{x}{\code{SpatRaster}. Preferably backed by a file so that workers can re-open it cheaply; in-memory rasters are shipped to the workers via \code{\link{wrap}}}
  \item{fun}{function. Takes a \code{SpatRaster} as its first argument and must return a \code{SpatRaster}. It is called once per tile, with the windowed \code{x} as the first argument and any extra \code{...} arguments appended}
  \item{cores}{one of: an integer (\code{> 1} creates a PSOCK cluster of that size, torn down at the end of the call); a \code{cluster} object created with \code{parallel::makeCluster}; the string \code{"future"} (use \code{\link[future.apply]{future_lapply}} with the currently-active \code{future::\link[future]{plan}}); or a future strategy object (e.g. \code{future::multisession()}), which is set as the plan for the duration of the call. The \code{"future"} and plan-object paths require the \pkg{future} and \pkg{future.apply} packages}
  \item{cpkgs}{character. Names of packages that must be loaded on each worker before \code{fun} is called (e.g. \code{"randomForest"})}
  \item{tiles}{specification of the tiles. If \code{NULL} (the default), the tile size is chosen automatically and a per-worker memory budget that depends on \code{cores} (see \code{\link{getTileExtents}}. Otherwise one of: a \code{SpatExtent}; a list of \code{SpatExtent}s; a 4-column matrix of extents \code{(xmin, xmax, ymin, ymax)} as returned by \code{\link{getTileExtents}}; a \code{SpatRaster} or \code{SpatVector} that defines the tile geometry; or one or two integers specifying the number of rows and columns \emph{per tile} (passed on to \code{\link{getTileExtents}})}
  \item{buffer}{integer. Number of additional rows and columns (one number, or two for separate row/column buffers) of \code{x} to read on each side of every tile, so that \code{fun} can compute neighbourhood-dependent values (e.g. \code{\link{focal}}) without edge artefacts. The result of \code{fun} is then cropped back to the un-buffered tile extent before it is written, so the per-tile outputs are non-overlapping. Only used when \code{tiles = NULL}; ignored with a warning otherwise (when supplying \code{tiles} explicitly, build the overlap into the tiles yourself with \code{getTileExtents(..., buffer=)} and pick an \code{overlap_fun})}
  \item{...}{additional arguments passed to \code{fun}. When \code{cores > 1} these are serialized and shipped to each worker; \code{SpatRaster}, \code{SpatVector} and other \pkg{terra} objects in \code{...} are automatically wrapped and unwrapped}
  \item{filename}{character. Output filename for the assembled result. When empty (the default) and \code{overlap_fun = NULL}, a virtual raster (VRT) is returned that references the per-tile files (see \emph{Details})}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten if it exists}
  \item{wopt}{list. Writing options as in \code{\link{writeRaster}}}
  \item{overlap_fun}{character or \code{NULL}. The function name used by \code{\link{mosaic}} to combine overlapping cells (e.g. \code{"mean"}, \code{"min"}, \code{"max"}, \code{"sum"}, \code{"first"}). Use this when \code{tiles} was built with overlap (e.g. \code{getTileExtents(..., buffer=)}). When \code{NULL} (the default), tiles are assembled with \code{\link{vrt}}, which is correct for non-overlapping tiles (including the auto path with \code{buffer}, which crops away the overlap) and is essentially free}
}

\details{
\strong{Buffered tiles for focal-style operations.} When \code{tiles = NULL} and \code{buffer > 0}, each tile is read on an expanded extent of (\code{buffer} cells on each side, clamped to \code{x}'s extent), \code{fun} is applied, and the result is cropped back to the un-buffered tile extent before it is written. This avoids the edge effects that operations like \code{\link{focal}} would otherwise produce at tile boundaries.

The per-tile files are then assembled into the final \code{SpatRaster}:
\itemize{
  \item When \code{overlap_fun = NULL} (the default), the assembly is a virtual raster (\code{\link{vrt}}). This is fast and correct as long as the tiles do not overlap. If they do, the value of the last tile in the overlapping regions is used.
  \item When \code{overlap_fun} is set (e.g. \code{"mean"}), the assembly uses \code{\link{mosaic}} and applies the function to the values of overlapping cells.
}

If \code{filename} is empty and a VRT is returned, the tiled data files are kept for the rest of the R session, but they are lost after the session.

Extra arguments passed via \code{...} must be named when \code{cores > 1}.

}

\value{
A \code{SpatRaster}.
}

\seealso{
\code{\link{getTileExtents}}, \code{\link{makeTiles}}, \code{\link{mosaic}}, \code{\link{vrt}}, \code{\link{window}}, \code{\link{wrap}}
}

\examples{
\donttest{
f <- system.file("ex/elev.tif", package="terra")
r <- rast(f)

# not parallel
# Returned object is a VRT backed by per-tile files in tempdir().
out1 <- tile_apply(r, function(x) x * 2)

# focal, use a buffer wide enough for the focal window
out3 <- tile_apply(r, function(x) focal(x, w=5, fun="mean", 
		na.rm=TRUE), buffer=2)

# parallel on 2 workers
\dontrun{
out4 <- tile_apply(r, function(x) x * 2, cores=2)

# pass extra arguments by name; nested terra objects are auto-wrapped
out5 <- tile_apply(r, function(x, k) x * k, cores=2, tiles=50, k=10)

# use a future plan instead of a cluster
if (requireNamespace("future", quietly=TRUE) &&
    requireNamespace("future.apply", quietly=TRUE)) {
    future::plan(future::multisession, workers=2)
    out6 <- tile_apply(r, function(x) x * 2, cores="future")
    future::plan(future::sequential)
}
}
}
}

\keyword{spatial}
