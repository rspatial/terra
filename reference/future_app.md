# Run app() or focal() over tiles using a future plan

Thin convenience wrappers around
[`tile_apply`](https://rspatial.github.io/terra/reference/tile_apply.md)
that route the per-tile work through the future framework via
[`future_lapply`](https://future.apply.futureverse.org/reference/future_lapply.html).
They:

- tile the input `SpatRaster` (auto-tiles by default),

- ship `x` to the workers with
  [`wrap`](https://rspatial.github.io/terra/reference/wrap.md) –
  file-backed sources travel as paths only; large in-memory rasters are
  materialised to a temporary file first, so very large rasters never
  have their values copied to every worker,

- run [`app`](https://rspatial.github.io/terra/reference/app.md)
  (`future_app`) or
  [`focal`](https://rspatial.github.io/terra/reference/focal.md)
  (`future_focal`) on each tile,

- assemble the per-tile outputs back into a single `SpatRaster`.

The user controls parallelism by setting up a
`future::`[`plan`](https://future.futureverse.org/reference/plan.html)
(e.g.
[`future::multisession`](https://future.futureverse.org/reference/multisession.html),
[`future::cluster`](https://future.futureverse.org/reference/cluster.html),
`future.batchtools`, ...). With the default
[`future::sequential`](https://future.futureverse.org/reference/sequential.html)
plan the call still works, but runs in-process.

## Usage

``` r
future_app(x, fun, ..., plan=NULL, cpkgs=NULL,
           tiles=NULL, buffer=0,
           filename="", overwrite=FALSE, wopt=list())

future_focal(x, w=3, fun="sum", ..., plan=NULL, cpkgs=NULL,
             tiles=NULL, buffer=NULL,
             filename="", overwrite=FALSE, wopt=list())
```

## Arguments

- x:

  `SpatRaster`

- fun:

  a function as accepted by
  [`app`](https://rspatial.github.io/terra/reference/app.md) (for
  `future_app`) or by
  [`focal`](https://rspatial.github.io/terra/reference/focal.md) (for
  `future_focal`). Forwarded verbatim to the per-tile call

- w:

  focal window size or matrix; passed to
  [`focal`](https://rspatial.github.io/terra/reference/focal.md)

- ...:

  additional arguments forwarded to the per-tile
  [`app()`](https://rspatial.github.io/terra/reference/app.md) /
  [`focal()`](https://rspatial.github.io/terra/reference/focal.md) call
  (e.g. `na.rm=`, `fillvalue=`, `cores=` for an inner per-tile
  parallelism, `wopt=`, ...)

- plan:

  optional `future` plan or strategy object (e.g.
  [`future::multisession()`](https://future.futureverse.org/reference/multisession.html)).
  When `NULL` (the default), the currently-active
  [`future::plan()`](https://future.futureverse.org/reference/plan.html)
  is used; the call does not modify it. When given, it is installed for
  the duration of the call and the previous plan is restored on exit

- cpkgs:

  character. Names of additional packages that workers must load (e.g.
  `"randomForest"`). Passed to `future.apply` as
  `future.packages = c("terra", cpkgs)`

- tiles:

  tile specification, see
  [`tile_apply`](https://rspatial.github.io/terra/reference/tile_apply.md).
  Defaults to `NULL` (auto tiles from GDAL block size + worker count)

- buffer:

  integer read buffer (in cells) used when `tiles = NULL`, see
  [`tile_apply`](https://rspatial.github.io/terra/reference/tile_apply.md).
  For `future_focal` the default is `floor(max(w)/2)`, which is wide
  enough that each per-tile
  [`focal()`](https://rspatial.github.io/terra/reference/focal.md) sees
  the same neighbours as a full-raster
  [`focal()`](https://rspatial.github.io/terra/reference/focal.md); for
  `future_app` the default of `0` is correct because
  [`app()`](https://rspatial.github.io/terra/reference/app.md) is
  per-cell

- filename:

  output filename; empty (default) returns a VRT-backed `SpatRaster`,
  see
  [`tile_apply`](https://rspatial.github.io/terra/reference/tile_apply.md)

- overwrite:

  logical

- wopt:

  list. Writing options as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.md)

## Details

These are equivalent to:

    tile_apply(x, function(z, ...) app(z, fun, ...),   ..., cores = "future")
    tile_apply(x, function(z, ...) focal(z, w, fun, ...), ...,
               buffer = floor(max(w)/2), cores = "future")

with
[`wrap()`](https://rspatial.github.io/terra/reference/wrap.md)/[`unwrap()`](https://rspatial.github.io/terra/reference/wrap.md)
of `x` and any terra objects in `...` done automatically. See
[`tile_apply`](https://rspatial.github.io/terra/reference/tile_apply.md)
for the full machinery (auto-tiling, buffered tiles, VRT vs.\\ `mosaic`
assembly, and the disk-streaming guarantee).

The wrappers are intentionally minimal so they compose with the rest of
the future ecosystem (`handlers`,
[`plan`](https://future.futureverse.org/reference/plan.html) stacks,
remote workers, ...).

## Value

A `SpatRaster`.

## See also

[`tile_apply`](https://rspatial.github.io/terra/reference/tile_apply.md),
[`app`](https://rspatial.github.io/terra/reference/app.md),
[`focal`](https://rspatial.github.io/terra/reference/focal.md),
[`wrap`](https://rspatial.github.io/terra/reference/wrap.md),
[`getTileExtents`](https://rspatial.github.io/terra/reference/makeTiles.md),
`future::`[`plan`](https://future.futureverse.org/reference/plan.html),
`future.apply::`[`future_lapply`](https://future.apply.futureverse.org/reference/future_lapply.html)

## Examples

``` r
if (FALSE) { # \dontrun{
if (requireNamespace("future", quietly=TRUE) &&
    requireNamespace("future.apply", quietly=TRUE)) {

    library(future)
    library(future.apply)
    plan(multisession, workers = 2)

    f <- system.file("ex/elev.tif", package="terra")
    r <- rast(f)

    # app(): per-cell. No buffer needed.
    out_a <- future_app(r, fun = function(v) v * 2)

    # focal(): auto-buffer = floor(w/2) keeps tile edges identical to a
    # full-raster focal().
    out_f <- future_focal(r, w = 5, fun = "mean", na.rm = TRUE)

    # explicit plan, restored on exit
    out_b <- future_app(r, fun = function(v) v + 1,
                        plan = sequential())

    plan(sequential)
}
} # }
```
