# Apply a function to tiles of a SpatRaster in parallel

A high-level convenience wrapper that splits a `SpatRaster` into tiles,
applies a user function to each tile (with
[`window`](https://rspatial.github.io/terra/reference/window.md) set
under the hood so no values are copied), and assembles the per-tile
results back into a single `SpatRaster`. The work can be distributed
over a cluster of worker processes.

This is the “Strategy 3” pattern from the terra parallelization guide,
packaged as one call.

## Usage

``` r
tile_apply(x, fun, cores=1, cpkgs=NULL, tiles=NULL, buffer=0, ...,
           filename="", overwrite=FALSE, wopt=list(),
           overlap_fun=NULL)
```

## Arguments

- x:

  `SpatRaster`. Should normally be backed by a file (so that workers can
  re-open it cheaply); in-memory rasters are shipped to the workers via
  [`wrap`](https://rspatial.github.io/terra/reference/wrap.md)

- fun:

  function. Takes a `SpatRaster` as its first argument and must return a
  `SpatRaster`. It is called once per tile, with the windowed `x` as the
  first argument and any extra `...` arguments appended

- cores:

  integer, or a `cluster` object created with
  [`parallel::makeCluster`](https://rdrr.io/r/parallel/makeCluster.html).
  If a number greater than 1, a PSOCK cluster of that size is created
  and torn down inside the call

- cpkgs:

  character. Names of packages that must be loaded on each worker before
  `fun` is called (e.g. `"randomForest"`)

- tiles:

  specification of the tiles. If `NULL` (the default), the tile size is
  chosen automatically from the GDAL block size of `x` and a per-worker
  memory budget that depends on `cores` (see
  [`getTileExtents`](https://rspatial.github.io/terra/reference/makeTiles.md)
  *Details*). Otherwise one of: a `SpatExtent`; a list of `SpatExtent`s;
  a 4-column matrix of extents `(xmin, xmax, ymin, ymax)` as returned by
  [`getTileExtents`](https://rspatial.github.io/terra/reference/makeTiles.md);
  a `SpatRaster` or `SpatVector` that defines the tile geometry; or one
  or two integers specifying the number of rows and columns *per tile*
  (passed on to
  [`getTileExtents`](https://rspatial.github.io/terra/reference/makeTiles.md))

- buffer:

  integer. Number of additional rows and columns (one number, or two for
  separate row/column buffers) of `x` to read on each side of every
  tile, so that `fun` can compute neighbourhood-dependent values (e.g.
  [`focal`](https://rspatial.github.io/terra/reference/focal.md))
  without edge artefacts. The result of `fun` is then cropped back to
  the un-buffered tile extent before it is written, so the per-tile
  outputs are non-overlapping and assemble cleanly. Only used when
  `tiles = NULL`; ignored with a warning otherwise (when supplying
  `tiles` explicitly, build the overlap into the tiles yourself with
  `getTileExtents(..., buffer=)` and pick an `overlap_fun`)

- ...:

  additional arguments passed to `fun`. When `cores > 1` these are
  serialized and shipped to each worker; `SpatRaster`, `SpatVector` and
  other terra objects in `...` are automatically wrapped and unwrapped

- filename:

  character. Output filename for the assembled result. When empty (the
  default) and `overlap_fun = NULL`, a virtual raster (VRT) is returned
  that references the per-tile files (see *Details*)

- overwrite:

  logical. If `TRUE`, `filename` is overwritten if it exists

- wopt:

  list. Writing options as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md),
  applied both to the per-tile intermediate files and to the final
  output

- overlap_fun:

  character or `NULL`. The function name used by
  [`mosaic`](https://rspatial.github.io/terra/reference/mosaic.md) to
  combine overlapping cells (e.g. `"mean"`, `"min"`, `"max"`, `"sum"`,
  `"first"`). Use this when `tiles` was built with overlap (e.g.
  `getTileExtents(..., buffer=)`). When `NULL` (the default), tiles are
  assembled with
  [`vrt`](https://rspatial.github.io/terra/reference/vrt.md), which is
  correct for non-overlapping tiles (including the auto path with
  `buffer`, which crops away the overlap) and is essentially free

## Details

For each tile, `tile_apply` sets a read window on `x`
(`window(x) <- e`), so no pixel values are copied: subsequent operations
on the windowed raster behave as if it had only that extent. The user
function `fun` therefore sees a normal `SpatRaster` restricted to the
tile.

**Memory safety.** Each tile's result is written to its own file in a
per-call sub-directory of
[`tempdir()`](https://rdrr.io/r/base/tempfile.html) (as a GeoTIFF,
honouring `wopt`). Workers therefore return only a filename, not pixel
data, so the parent process never holds more than one tile's worth of
values in RAM. The sequential path follows the same contract.

**Buffered tiles for focal-style operations.** When `tiles = NULL` and
`buffer > 0`, each tile is read on a slightly expanded extent (`buffer`
cells on each side, clamped to `x`'s extent), `fun` is applied, and the
result is cropped back to the un-buffered tile extent before it is
written. This avoids the edge effects that operations like
[`focal`](https://rspatial.github.io/terra/reference/focal.md) would
otherwise produce at tile boundaries, while still letting the assembled
output use the cheap, lossless
[`vrt`](https://rspatial.github.io/terra/reference/vrt.md) path. The
buffer also enlarges the per-worker peak RAM (a tile reads
`(rows + 2 buffer)(cols + 2 buffer)` cells), which is usually negligible
since `buffer` is small compared with the auto tile dimensions, but if
you choose a large buffer you may want to raise `cores` so the auto
sizer picks smaller tiles.

The per-tile files are then assembled into the final `SpatRaster`:

- When `overlap_fun = NULL` (the default), the assembly is a virtual
  raster ([`vrt`](https://rspatial.github.io/terra/reference/vrt.md)).
  This is essentially free and never reads the pixel data again. It is
  correct as long as the tiles do not overlap; if they do (e.g. when
  built with `getTileExtents(..., buffer=)`), GDAL keeps the value of
  the last tile drawn in the overlap region.

- When `overlap_fun` is set (e.g. `"mean"`), the assembly uses
  [`mosaic`](https://rspatial.github.io/terra/reference/mosaic.md) which
  streams the data through terra's standard chunked I/O and applies
  `overlap_fun` on overlapping cells. This is the right path for
  buffered tiles, focal operations, and similar.

If `filename` is supplied, the assembled result is materialised there
and the per-tile intermediate files are removed. If `filename` is empty
and a VRT is being built, the per-tile files are kept for the rest of
the R session because they back the returned raster.

When `cores > 1`, `x` is shipped to the workers with
[`wrap`](https://rspatial.github.io/terra/reference/wrap.md). For
file-backed rasters this is cheap (only the metadata and filename travel
across the worker boundary, and each worker re-opens the file
independently). For in-memory rasters the values are copied to each
worker.

Extra arguments passed via `...` must be named when `cores > 1`.

## Value

A `SpatRaster`.

## See also

[`getTileExtents`](https://rspatial.github.io/terra/reference/makeTiles.md),
[`makeTiles`](https://rspatial.github.io/terra/reference/makeTiles.md),
[`mosaic`](https://rspatial.github.io/terra/reference/mosaic.md),
[`vrt`](https://rspatial.github.io/terra/reference/vrt.md),
[`window`](https://rspatial.github.io/terra/reference/window.md),
[`wrap`](https://rspatial.github.io/terra/reference/wrap.md),
[`app`](https://rspatial.github.io/terra/reference/app.md)

## Examples

``` r
# \donttest{
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
if (FALSE) { # \dontrun{
out4 <- tile_apply(r, function(x) x * 2, cores=2)

# pass extra arguments by name; nested terra objects are auto-wrapped
out5 <- tile_apply(r, function(x, k) x * k, cores=2, tiles=50, k=10)
} # }
# }
```
