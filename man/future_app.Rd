\name{future_app}

\alias{future_app}
\alias{future_focal}

\title{Run app() or focal() over tiles using a future plan}

\description{
Thin convenience wrappers around \code{\link{tile_apply}} that route the
per-tile work through the \pkg{future} framework via
\code{\link[future.apply]{future_lapply}}. They:

\itemize{
  \item tile the input \code{SpatRaster} (auto-tiles by default),
  \item ship \code{x} to the workers with \code{\link{wrap}} -- file-backed
        sources travel as paths only; large in-memory rasters are
        materialised to a temporary file first, so very large rasters
        never have their values copied to every worker,
  \item run \code{\link{app}} (\code{future_app}) or \code{\link{focal}}
        (\code{future_focal}) on each tile,
  \item assemble the per-tile outputs back into a single \code{SpatRaster}.
}

The user controls parallelism by setting up a \code{future::\link[future]{plan}}
(e.g. \code{future::multisession}, \code{future::cluster},
\code{future.batchtools}, ...). With the default \code{future::sequential}
plan the call still works, but runs in-process.
}

\usage{
future_app(x, fun, ..., plan=NULL, cpkgs=NULL,
           tiles=NULL, buffer=0,
           filename="", overwrite=FALSE, wopt=list())

future_focal(x, w=3, fun="sum", ..., plan=NULL, cpkgs=NULL,
             tiles=NULL, buffer=NULL,
             filename="", overwrite=FALSE, wopt=list())
}

\arguments{
  \item{x}{\code{SpatRaster}}
  \item{fun}{a function as accepted by \code{\link{app}} (for \code{future_app}) or by \code{\link{focal}} (for \code{future_focal}). Forwarded verbatim to the per-tile call}
  \item{w}{focal window size or matrix; passed to \code{\link{focal}}}
  \item{...}{additional arguments forwarded to the per-tile \code{app()} / \code{focal()} call (e.g. \code{na.rm=}, \code{fillvalue=}, \code{cores=} for an inner per-tile parallelism, \code{wopt=}, ...)}
  \item{plan}{optional \code{future} plan or strategy object (e.g. \code{future::multisession()}). When \code{NULL} (the default), the currently-active \code{future::plan()} is used; the call does not modify it. When given, it is installed for the duration of the call and the previous plan is restored on exit}
  \item{cpkgs}{character. Names of additional packages that workers must load (e.g. \code{"randomForest"}). Passed to \code{future.apply} as \code{future.packages = c("terra", cpkgs)}}
  \item{tiles}{tile specification, see \code{\link{tile_apply}}. Defaults to \code{NULL} (auto tiles from GDAL block size + worker count)}
  \item{buffer}{integer read buffer (in cells) used when \code{tiles = NULL}, see \code{\link{tile_apply}}. For \code{future_focal} the default is \code{floor(max(w)/2)}, which is wide enough that each per-tile \code{focal()} sees the same neighbours as a full-raster \code{focal()}; for \code{future_app} the default of \code{0} is correct because \code{app()} is per-cell}
  \item{filename}{output filename; empty (default) returns a VRT-backed \code{SpatRaster}, see \code{\link{tile_apply}}}
  \item{overwrite}{logical}
  \item{wopt}{list. Writing options as in \code{\link{writeRaster}}}
}

\details{
These are equivalent to:

\preformatted{
tile_apply(x, function(z, ...) app(z, fun, ...),   ..., cores = "future")
tile_apply(x, function(z, ...) focal(z, w, fun, ...), ...,
           buffer = floor(max(w)/2), cores = "future")
}

with \code{wrap()}/\code{unwrap()} of \code{x} and any \pkg{terra} objects
in \code{...} done automatically. See \code{\link{tile_apply}} for the full
machinery (auto-tiling, buffered tiles, VRT vs.\ \code{mosaic} assembly,
and the disk-streaming guarantee).

The wrappers are intentionally minimal so they compose with the rest of
the \pkg{future} ecosystem (\code{\link[progressr]{handlers}},
\code{\link[future]{plan}} stacks, remote workers, ...).
}

\value{A \code{SpatRaster}.}

\seealso{
\code{\link{tile_apply}}, \code{\link{app}}, \code{\link{focal}},
\code{\link{wrap}}, \code{\link{getTileExtents}},
\code{future::\link[future]{plan}},
\code{future.apply::\link[future.apply]{future_lapply}}
}

\examples{
\dontrun{
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
}
}

\keyword{spatial}
