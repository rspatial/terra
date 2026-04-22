# tests for future integration: tile_apply(cores="future"|plan), future_app(),
# future_focal(). All gated on the optional packages being installed so the
# CRAN check passes when they are not.

if (!requireNamespace("future", quietly=TRUE) ||
	!requireNamespace("future.apply", quietly=TRUE)) {
	exit_file("future / future.apply not installed; skipping")
}

f <- system.file("ex/elev.tif", package="terra")
r <- rast(f)
v_in <- values(r)


# helper: install a sequential plan for the duration of the file (so we
# never spawn worker processes and stay well under R CMD check --as-cran's
# 2-core limit). Restored on exit().
old_plan <- future::plan(future::sequential)
on.exit(future::plan(old_plan), add=TRUE)


# tile_apply(cores="future") with the active sequential plan: equivalent
# to a single-process tile_apply() but exercising the future.apply path.

out_seq <- tile_apply(r, function(x) x * 2, tiles=50, cores=1L)
out_fut <- tile_apply(r, function(x) x * 2, tiles=50, cores="future")
expect_inherits(out_fut, "SpatRaster")
expect_equal(dim(out_fut), dim(r))
v_seq <- values(out_seq); v_fut <- values(out_fut)
ok <- !is.na(v_seq) & !is.na(v_fut)
expect_true(sum(ok) > 0)
expect_equal(v_seq[ok], v_fut[ok])


# tile_apply with a future strategy passed via cores= (no parens).
# This installs the plan for the duration of the call and restores it.
out_plan <- tile_apply(r, function(x) x + 1, tiles=50,
					   cores=future::sequential)
expect_inherits(out_plan, "SpatRaster")
v <- values(out_plan)
ok <- !is.na(v) & !is.na(v_in)
expect_equal(v[ok], v_in[ok] + 1)
# the plan we set up-front should be restored after the call
expect_inherits(future::plan(), "sequential")


# extra `...` args travel across the future boundary (named, like cluster path)
out_dots <- tile_apply(r, function(x, k) x * k, tiles=50, cores="future", k=4)
v <- values(out_dots)
ok <- !is.na(v) & !is.na(v_in)
expect_equal(v[ok], v_in[ok] * 4)


# future_app(): per-cell, no buffer needed, must equal app() on the full raster
out_app <- future_app(r, fun = function(v) v * 2 + 1)
ref_app <- app(r,    fun = function(v) v * 2 + 1)
v_app <- values(out_app); v_ref <- values(ref_app)
ok <- !is.na(v_app) & !is.na(v_ref)
expect_equal(v_app[ok], v_ref[ok])


# future_focal() with the auto buffer matches a full-raster focal exactly
out_fc <- future_focal(r, w=5, fun="mean", na.rm=TRUE,
					   wopt=list(datatype="FLT8S"))
ref_fc <- focal(r, w=5, fun="mean", na.rm=TRUE)
v_fc <- values(out_fc); v_ref <- values(ref_fc)
ok <- !is.na(v_fc) & !is.na(v_ref)
expect_true(sum(ok) > 0)
expect_equal(v_fc[ok], v_ref[ok])


# future_focal() with explicit tiles disables the auto buffer cleanly
tpl <- rast(ext(r), nrows=2, ncols=2, crs=crs(r))
out_fc2 <- future_focal(r, w=3, fun="mean", na.rm=TRUE, tiles=tpl,
						wopt=list(datatype="FLT8S"))
expect_inherits(out_fc2, "SpatRaster")
expect_equal(dim(out_fc2), dim(r))


# future_app() to a filename
ftif <- file.path(tempdir(), "future_app_test.tif")
unlink(ftif)
out_file <- future_app(r, fun = function(v) v + 10, tiles=50,
					   filename=ftif, overwrite=TRUE,
					   wopt=list(datatype="FLT8S"))
expect_true(file.exists(ftif))
v <- values(out_file)
ok <- !is.na(v) & !is.na(v_in)
expect_equal(v[ok], v_in[ok] + 10)
unlink(ftif)


# future_app() with an explicit plan= argument (strategy without parens):
# the plan is installed only for the duration of the call.
out_plan2 <- future_app(r, fun = function(v) v - 5,
						tiles=50, plan=future::sequential)
v <- values(out_plan2)
ok <- !is.na(v) & !is.na(v_in)
expect_equal(v[ok], v_in[ok] - 5)
expect_inherits(future::plan(), "sequential")
