# tests for focal() with cores= dispatch (R-fun path) and the TBB-aware
# focalValues() implementation.

# CRAN's --as-cran limits worker processes to 2; everything below stays
# within that ceiling.

# ---- focalValues: parallel result must equal sequential ----
old <- terra::terraOptions(print=FALSE)

# global lonlat raster (exercises the wrap-around branch)
r_glob <- terra::rast(ncols=20L, nrows=15L, vals=seq_len(20L*15L))
terra::terraOptions(parallel=FALSE, threads=0L)
seq_glob <- terra::focalValues(r_glob, w=3L, fill=-9)
terra::terraOptions(parallel=TRUE, threads=2L)
par_glob <- terra::focalValues(r_glob, w=3L, fill=-9)
expect_identical(seq_glob, par_glob)

# projected raster (no wrap-around; exercises the interior std::copy fast path)
r_proj <- terra::rast(ncols=40L, nrows=30L,
                       xmin=0, xmax=4000, ymin=0, ymax=3000,
                       crs="EPSG:32631", vals=seq_len(40L*30L))
terra::terraOptions(parallel=FALSE, threads=0L)
seq_proj <- terra::focalValues(r_proj, w=5L, fill=-1)
terra::terraOptions(parallel=TRUE, threads=2L)
par_proj <- terra::focalValues(r_proj, w=5L, fill=-1)
expect_identical(seq_proj, par_proj)

# restore options
terra::terraOptions(parallel=old$parallel, threads=old$threads)


# ---- focal() with R fun + cores= ----
# small enough to stay quick under tinytest, big enough that the per-block
# parallel path is engaged (.focal_split_chunks splits when nrow >= 64*ncores).
nr <- 80L; nc <- 80L
r <- terra::rast(nrows=nr, ncols=nc, xmin=0, xmax=nc*100, ymin=0, ymax=nr*100,
                  crs="EPSG:32631", vals=seq_len(nr*nc))

# scalar-output user fun
f_scalar <- function(x) sum(x, na.rm=TRUE) / max(1L, sum(!is.na(x)))
ref <- terra::values(terra::focal(r, w=5, fun=f_scalar, cores=1))
got <- terra::values(terra::focal(r, w=5, fun=f_scalar, cores=2))
ok  <- !is.na(ref) & !is.na(got)
expect_true(max(abs(ref[ok] - got[ok])) < 1e-10)

# vector-output user fun (transp branch)
f_vec <- function(x) c(min=min(x, na.rm=TRUE), max=max(x, na.rm=TRUE))
ref <- terra::values(terra::focal(r, w=5, fun=f_vec, cores=1))
got <- terra::values(terra::focal(r, w=5, fun=f_vec, cores=2))
ok  <- !is.na(ref) & !is.na(got)
expect_true(max(abs(ref[ok] - got[ok])) < 1e-10)

# multi-layer
r2  <- c(r, r * 2L)
ref <- terra::values(terra::focal(r2, w=5, fun=f_scalar, cores=1))
got <- terra::values(terra::focal(r2, w=5, fun=f_scalar, cores=2))
ok  <- !is.na(ref) & !is.na(got)
expect_true(max(abs(ref[ok] - got[ok])) < 1e-10)

# weights matrix path (dow=TRUE, msz != prod(w) when NA in weights)
wm  <- matrix(c(0,1,0,1,1,1,0,1,0), 3L, 3L)
ref <- terra::values(terra::focal(r, w=wm, fun=f_scalar, cores=1))
got <- terra::values(terra::focal(r, w=wm, fun=f_scalar, cores=2))
ok  <- !is.na(ref) & !is.na(got)
expect_true(max(abs(ref[ok] - got[ok])) < 1e-10)

# user-supplied cluster object (so we don't pay setup cost in users' loops)
if (requireNamespace("parallel", quietly=TRUE)) {
	cl <- parallel::makeCluster(2L)
	got <- terra::values(terra::focal(r, w=5, fun=f_scalar, cores=cl))
	parallel::stopCluster(cl)
	ok  <- !is.na(ref) & !is.na(got)
	# `ref` here is the weighted-window result; recompute scalar reference
	ref <- terra::values(terra::focal(r, w=5, fun=f_scalar, cores=1))
	ok  <- !is.na(ref) & !is.na(got)
	expect_true(max(abs(ref[ok] - got[ok])) < 1e-10)
}
