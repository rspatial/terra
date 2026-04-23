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

