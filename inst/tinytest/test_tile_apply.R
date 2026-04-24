# tests for tile_apply()

f <- system.file("ex/elev.tif", package="terra")
r <- rast(f)


# numeric tile spec: x * 2 should reproduce r * 2 cell-for-cell 

out2 <- tile_apply(r, function(x) x * 2, tiles=50)
expect_inherits(out2, "SpatRaster")
expect_equal(dim(out2), dim(r))

v_in <- values(r)
v_out <- values(out2)
ok <- !is.na(v_in) & !is.na(v_out)
expect_true(sum(ok) > 0)
expect_equal(v_out[ok], 2 * v_in[ok])


# tiles given as a 4-column matrix (from getTileExtents) 

m <- getTileExtents(r, 50)
out_m <- tile_apply(r, function(x) x + 1, tiles=m)
expect_equal(dim(out_m), dim(r))
v <- values(out_m)
ok <- !is.na(v) & !is.na(v_in)
expect_equal(v[ok], v_in[ok] + 1)


# tiles given as a list of SpatExtents

elist <- lapply(seq_len(nrow(m)), function(i) ext(m[i, ]))
out_list <- tile_apply(r, function(x) x + 1, tiles=elist)
expect_equal(dim(out_list), dim(r))
v <- values(out_list)
ok <- !is.na(v) & !is.na(v_in)
expect_equal(v[ok], v_in[ok] + 1)


# tiles given as a SpatRaster template

tpl <- rast(ext(r), nrows=2, ncols=2, crs=crs(r))
out_tpl <- tile_apply(r, function(x) x * 2, tiles=tpl)
expect_equal(dim(out_tpl), dim(r))


# extra `...` arguments are forwarded by name

out_dots <- tile_apply(r, function(x, k) x * k, tiles=50, k=3)
v <- values(out_dots)
ok <- !is.na(v) & !is.na(v_in)
expect_equal(v[ok], v_in[ok] * 3)


# single-tile shortcut (no mosaic call) 

out_one <- tile_apply(r, function(x) x + 0, tiles=ext(r))
expect_equal(dim(out_one), dim(r))


# focal: buffered tiles + overlap_fun = "mean" assembles via mosaic ─
# (per-tile output includes buffered cells; mean blending across overlap.)

tiles_buf <- getTileExtents(r, 50, buffer=5)
out_focal <- tile_apply(r,
	function(x) focal(x, w=5, fun="mean", na.rm=TRUE),
	tiles=tiles_buf, overlap_fun="mean")
expect_equal(dim(out_focal), dim(r))
expect_true(any(!is.na(values(out_focal))))


# memory-safe assembly: default returns a VRT-backed SpatRaster 
out_vrt <- tile_apply(r, function(x) x * 2, tiles=50)
expect_equal(as.vector(minmax(out_vrt)), c(282, 1094))

v <- values(out_vrt)
ok <- !is.na(v) & !is.na(v_in)
expect_equal(as.vector(table(ok)), c(3942, 4608))

expect_true(grepl("\\.vrt$", sources(out_vrt)[1]))


# memory-safe assembly: filename materialises a single output file

ftif_seq <- file.path(tempdir(), "tile_apply_test_seq.tif")
unlink(ftif_seq)
out_file <- tile_apply(r, function(x) x + 1, tiles=50,
	filename=ftif_seq, overwrite=TRUE)
expect_true(file.exists(ftif_seq))
expect_equal(normalizePath(sources(out_file)[1], "/"),
			 normalizePath(ftif_seq, "/"))
v <- values(out_file)
ok <- !is.na(v) & !is.na(v_in)
expect_equal(v[ok], v_in[ok] + 1)
unlink(ftif_seq)


# overlap_fun = "first" via mosaic equals input on non-overlapping tiles

out_first <- tile_apply(r, function(x) x + 0, tiles=50, overlap_fun="first")
expect_equal(dim(out_first), dim(r))
v <- values(out_first)
ok <- !is.na(v) & !is.na(v_in)
expect_equal(v[ok], v_in[ok])


# auto buffer: focal(buffer=2) matches a full-raster focal exactly ──
# (per-tile output is cropped to the un-buffered tile extent, so the VRT path
#  yields the same answer as running focal on the whole raster. We write the
#  intermediates as FLT8S so the comparison is bit-exact; with the default
#  FLT4S a ~1e-5 roundtrip difference is normal.)
# cores=2 (not 4) so we stay within R CMD check --as-cran's
# _R_CHECK_LIMIT_CORES_ ceiling of 2 worker processes.

ref_focal <- focal(r, w=5, fun="mean", na.rm=TRUE)
out_buf <- tile_apply(r,
	function(x) focal(x, w=5, fun="mean", na.rm=TRUE),
	buffer=2, cores=2, wopt=list(datatype="FLT8S"))
expect_equal(dim(out_buf), dim(r))
v_ref <- values(ref_focal); v_out <- values(out_buf)
ok <- !is.na(v_ref) & !is.na(v_out)
expect_equal(v_out[ok], v_ref[ok])

# without buffer, multi-tile auto path produces seam artefacts at tile edges
out_nobuf <- tile_apply(r,
	function(x) focal(x, w=5, fun="mean", na.rm=TRUE),
	buffer=0, cores=2, wopt=list(datatype="FLT8S"))
v_nb <- values(out_nobuf)
ok <- !is.na(v_ref) & !is.na(v_nb)
expect_true(max(abs(v_ref[ok] - v_nb[ok])) > 0)


# auto buffer: works in parallel and via filename

ftif_buf <- file.path(tempdir(), "tile_apply_test_buf.tif")
unlink(ftif_buf)
out_bp <- tile_apply(r,
	function(x) focal(x, w=5, fun="mean", na.rm=TRUE),
	buffer=2, cores=2, filename=ftif_buf, overwrite=TRUE,
	wopt=list(datatype="FLT8S"))
expect_true(file.exists(ftif_buf))
v_bp <- values(out_bp); ok <- !is.na(v_ref) & !is.na(v_bp)
expect_equal(v_bp[ok], v_ref[ok])
unlink(ftif_buf)


# buffer is ignored (with a warning) when tiles is supplied explicitly

expect_warning(
	out_bw <- tile_apply(r, function(x) x + 1, tiles=50, buffer=3),
	pattern="buffer.*only used when 'tiles' is NULL"
)
v <- values(out_bw)
ok <- !is.na(v) & !is.na(v_in)
expect_equal(v[ok], v_in[ok] + 1)


# auto tiles: tiles=NULL uses GDAL block size + cores-aware sizing

# elev.tif's GDAL block is one strip (43 rows x 95 cols), so the auto path
# produces ~2 tiles regardless of cores (the one-block floor wins).
out_auto <- tile_apply(r, function(x) x * 2)
expect_inherits(out_auto, "SpatRaster")
expect_equal(dim(out_auto), dim(r))
v <- values(out_auto)
ok <- !is.na(v) & !is.na(v_in)
expect_equal(v[ok], 2 * v_in[ok])

# auto tiles with `cores>1` should still produce a correct mosaic
out_auto2 <- tile_apply(r, function(x) x + 1, cores=2)
expect_equal(dim(out_auto2), dim(r))
v <- values(out_auto2)
ok <- !is.na(v) & !is.na(v_in)
expect_equal(v[ok], v_in[ok] + 1)

# auto tiles: scaling with cores on a larger in-memory raster

big <- rast(nrows=2000, ncols=2000, vals=1)
m1 <- getTileExtents(big, cores=1)
m4 <- getTileExtents(big, cores=4)
expect_inherits(m1, "matrix")
expect_equal(ncol(m1), 4L)
# more workers must not produce fewer tiles
expect_true(nrow(m4) >= nrow(m1))
# at least cores * 1 tiles for parallelism
expect_true(nrow(m4) >= 4)

# auto tiles: align to GDAL block size on a tiled GeoTIFF

ftif <- file.path(tempdir(), "tile_apply_test_tiled.tif")
writeRaster(big, ftif, overwrite=TRUE,
	gdal=c("TILED=YES", "BLOCKXSIZE=128", "BLOCKYSIZE=128"))
if (file.exists(ftif)) {
	tiled <- rast(ftif)

	# fileBlocksize should report the tile size we set
	fb <- fileBlocksize(tiled)
	expect_equal(as.integer(fb[1, "rows"]), 128L)
	expect_equal(as.integer(fb[1, "cols"]), 128L)

	# Auto tile extents should align to whole 128x128 blocks. Interior (non-
	# edge) tiles all share the same nominal width/height; that nominal size
	# must be a multiple of 128. Edge tiles may be smaller.
	m_auto <- getTileExtents(tiled, cores=4)
	xres_t <- xres(tiled); yres_t <- yres(tiled)
	tile_w_cells <- round((m_auto[, "xmax"] - m_auto[, "xmin"]) / xres_t)
	tile_h_cells <- round((m_auto[, "ymax"] - m_auto[, "ymin"]) / yres_t)
	nominal_w <- max(tile_w_cells)
	nominal_h <- max(tile_h_cells)
	expect_equal(nominal_w %% 128L, 0L)
	expect_equal(nominal_h %% 128L, 0L)
	# and we should get at least as many tiles for more cores
	m_auto1 <- getTileExtents(tiled, cores=1)
	expect_true(nrow(m_auto) >= nrow(m_auto1))

	rm(tiled)
	unlink(ftif)
}


# error handling

expect_error(tile_apply("not a raster", function(x) x, tiles=50))
expect_error(tile_apply(r, function(x) values(x), tiles=50))   # fun must return a SpatRaster

# x with a window already set is rejected (use a window inside r's extent)
rw <- r
window(rw) <- ext(5.8, 6.0, 49.5, 49.8)
expect_true(window(rw))
expect_error(tile_apply(rw, function(x) x, tiles=50))


# parallel run ─
# Wrapped in tryCatch so a worker that cannot start does not break the suite.

if (requireNamespace("parallel", quietly=TRUE)) {
	tried <- tryCatch({
		out_par <- tile_apply(r, function(x) x * 2, tiles=50, cores=2)
		TRUE
	}, error = function(e) FALSE)
	if (tried) {
		v_par <- values(out_par)
		ok <- !is.na(v_par) & !is.na(v_in)
		expect_equal(v_par[ok], 2 * v_in[ok])

		# named extra args forwarded across workers
		out_par2 <- tile_apply(r, function(x, k) x * k, tiles=50, cores=2, k=7)
		v <- values(out_par2)
		ok <- !is.na(v) & !is.na(v_in)
		expect_equal(v[ok], 7 * v_in[ok])

		# user-supplied cluster object honored (not stopped by tile_apply)
		cl <- parallel::makeCluster(2)
		out_par3 <- tile_apply(r, function(x) x + 1, tiles=50, cores=cl)
		expect_inherits(out_par3, "SpatRaster")
		# cluster should still be alive after the call
		expect_silent(parallel::clusterCall(cl, function() 1L))
		parallel::stopCluster(cl)

		# parallel + filename: workers stream pixels straight to disk
		ftif_par <- file.path(tempdir(), "tile_apply_test_par.tif")
		unlink(ftif_par)
		out_par_f <- tile_apply(r, function(x) x * 3, tiles=50, cores=2,
			filename=ftif_par, overwrite=TRUE)
		expect_true(file.exists(ftif_par))
		v <- values(out_par_f)
		ok <- !is.na(v) & !is.na(v_in)
		expect_equal(v[ok], 3 * v_in[ok])
		unlink(ftif_par)


		# in-memory dispatch: when sources(x) are all empty, tile_apply
		# ships per-tile raster slices to workers instead of wrapping the
		# whole raster. Result must still match a sequential run cell-for-
		# cell, both with and without buffered tiles.
		mem_r <- rast(nrows=80, ncols=120, vals=runif(80*120))
		expect_true(all(sources(mem_r) == ""))

		out_mem_seq <- tile_apply(mem_r, function(x) x + 1, tiles=30,
			wopt=list(datatype="FLT8S"))
		out_mem_par <- tile_apply(mem_r, function(x) x + 1, tiles=30, cores=2,
			wopt=list(datatype="FLT8S"))
		expect_equal(as.vector(values(out_mem_par)),
					 as.vector(values(out_mem_seq)))

		# named extra args still forwarded on the slice path
		out_mem_par_k <- tile_apply(mem_r, function(x, k) x * k,
			tiles=30, cores=2, k=4, wopt=list(datatype="FLT8S"))
		expect_equal(as.vector(values(out_mem_par_k)),
					 as.vector(values(mem_r)) * 4)

		# buffered focal in parallel on an in-memory raster matches the
		# whole-raster focal exactly (FLT8S to avoid roundtrip noise).
		ref_focal_mem <- focal(mem_r, w=3, fun="mean", na.rm=TRUE)
		out_mem_buf <- tile_apply(mem_r,
			function(x) focal(x, w=3, fun="mean", na.rm=TRUE),
			buffer=1, cores=2, wopt=list(datatype="FLT8S"))
		v_ref <- values(ref_focal_mem); v_out <- values(out_mem_buf)
		ok <- !is.na(v_ref) & !is.na(v_out)
		expect_equal(v_out[ok], v_ref[ok])
	}
}
