
x <- y <- seq(0, 3, 1)
d <- terra::distance(cbind(x, y), lonlat=FALSE, sequential=TRUE)
sq2 <- sqrt(2)

expect_equal(d, c(0, rep(sq2, 3)))

d <- terra::distance(cbind(x, y), lonlat=FALSE, sequential=FALSE)
expect_equal(as.vector(d), c(sq2, 2*sq2, 3*sq2, sq2, 2*sq2, sq2))

d <- terra::distance(cbind(x, y), lonlat=TRUE, sequential=TRUE)
expect_equivalent(d, c(0.0, 156899.6, 156876.1, 156829.3), .000001)

d <- terra::distance(cbind(x, y), lonlat=TRUE, sequential=FALSE)
expect_equivalent(as.vector(d), c(156899.6, 313775.7, 470605.0, 156876.1, 313705.4, 156829.3), .000001)


# distance(SpatRaster, target=...) lon/lat fast path: must agree with the
# original implementation (gated by TERRA_DIST_SLOW=1) to within FP precision.
dist_both <- function(...) {
	old <- Sys.getenv("TERRA_DIST_SLOW", unset = NA)
	on.exit({
		if (is.na(old)) Sys.unsetenv("TERRA_DIST_SLOW")
		else Sys.setenv(TERRA_DIST_SLOW = old)
	})
	Sys.setenv(TERRA_DIST_SLOW = "0"); rf <- terra::distance(...)
	Sys.setenv(TERRA_DIST_SLOW = "1"); rs <- terra::distance(...)
	list(fast = terra::values(rf), slow = terra::values(rs))
}

set.seed(1)
rr <- terra::rast(nrows = 30, ncols = 60, xmin = -180, xmax = 180,
                  ymin = -90, ymax = 90, crs = "EPSG:4326")
terra::values(rr) <- 0
rr[c(150, 600, 1200, 1500)] <- 1

for (mth in c("haversine", "cosine")) {
	r <- dist_both(rr, target = 1, method = mth)
	expect_equal(is.na(r$fast), is.na(r$slow))
	fin <- !is.na(r$fast) & !is.na(r$slow)
	expect_true(all(abs(r$fast[fin] - r$slow[fin]) < 1e-4))   # < 0.1 mm
}

# target = NA path
rr2 <- terra::rast(nrows = 25, ncols = 50, xmin = -120, xmax = 120,
                   ymin = -60, ymax = 60, crs = "EPSG:4326")
terra::values(rr2) <- NA
rr2[c(40, 200, 800)] <- 5
r <- dist_both(rr2, method = "haversine")
expect_equal(is.na(r$fast), is.na(r$slow))
fin <- !is.na(r$fast) & !is.na(r$slow)
expect_true(all(abs(r$fast[fin] - r$slow[fin]) < 1e-4))


# maxdist: the fast path prunes in-loop using the threshold. It must produce
# the same NA pattern and the same finite values as the slow path (which
# filters after the fact).
for (md in c(5e6, 1e6, 2e5)) {
	r <- dist_both(rr, target = 1, method = "haversine", maxdist = md)
	expect_equal(is.na(r$fast), is.na(r$slow))
	fin <- !is.na(r$fast) & !is.na(r$slow)
	expect_true(all(abs(r$fast[fin] - r$slow[fin]) < 1e-4))
}

# maxdist in km
r <- dist_both(rr, target = 1, method = "haversine", unit = "km", maxdist = 500)
expect_equal(is.na(r$fast), is.na(r$slow))
fin <- !is.na(r$fast) & !is.na(r$slow)
expect_true(all(abs(r$fast[fin] - r$slow[fin]) < 1e-7))


# geo fast path: at EPSG:4326 both paths use WGS84, so they must be
# bit-equivalent (up to FP noise).
r <- dist_both(rr, target = 1, method = "geo")
expect_equal(is.na(r$fast), is.na(r$slow))
fin <- !is.na(r$fast) & !is.na(r$slow)
expect_true(all(abs(r$fast[fin] - r$slow[fin]) < 1e-6))

# geo + maxdist in meters
for (md in c(5e6, 1e6, 2e5)) {
	r <- dist_both(rr, target = 1, method = "geo", maxdist = md)
	expect_equal(is.na(r$fast), is.na(r$slow))
	fin <- !is.na(r$fast) & !is.na(r$slow)
	expect_true(all(abs(r$fast[fin] - r$slow[fin]) < 1e-6))
}

# geo + maxdist in km
r <- dist_both(rr, target = 1, method = "geo", unit = "km", maxdist = 500)
expect_equal(is.na(r$fast), is.na(r$slow))
fin <- !is.na(r$fast) & !is.na(r$slow)
expect_true(all(abs(r$fast[fin] - r$slow[fin]) < 1e-9))

# geo with target = NA (edges of non-NA region)
rr3 <- terra::rast(nrows = 25, ncols = 50, xmin = -120, xmax = 120,
                   ymin = -60, ymax = 60, crs = "EPSG:4326")
terra::values(rr3) <- NA
rr3[c(40, 200, 800)] <- 5
r <- dist_both(rr3, method = "geo")
expect_equal(is.na(r$fast), is.na(r$slow))
fin <- !is.na(r$fast) & !is.na(r$slow)
expect_true(all(abs(r$fast[fin] - r$slow[fin]) < 1e-6))


# distance(SpatRaster, SpatVector) lon/lat fast path for polygons/lines.
dist_both_v <- function(...) {
	old <- Sys.getenv("TERRA_DIST_VECTOR_SLOW", unset = NA)
	on.exit({
		if (is.na(old)) Sys.unsetenv("TERRA_DIST_VECTOR_SLOW")
		else Sys.setenv(TERRA_DIST_VECTOR_SLOW = old)
	})
	Sys.setenv(TERRA_DIST_VECTOR_SLOW = "0"); rf <- terra::distance(...)
	Sys.setenv(TERRA_DIST_VECTOR_SLOW = "1"); rs <- terra::distance(...)
	list(fast = terra::values(rf), slow = terra::values(rs))
}

rv <- terra::rast(nrows = 30, ncols = 60, xmin = -180, xmax = 180,
                  ymin = -60, ymax = 60, crs = "EPSG:4326")
poly <- terra::vect(matrix(c(-60, -10,  60, -10,  60, 40,  -60, 40,  -60, -10),
                           ncol = 2, byrow = TRUE),
                    type = "polygons", crs = "EPSG:4326")
line <- terra::vect(matrix(c(-120, -30,  -60, 0,  0, 30,  60, 0,  120, -30),
                           ncol = 2, byrow = TRUE),
                    type = "lines", crs = "EPSG:4326")

for (mth in c("haversine", "cosine", "geo")) {
	r <- dist_both_v(rv, poly, method = mth)
	expect_equal(is.na(r$fast), is.na(r$slow))
	fin <- !is.na(r$fast) & !is.na(r$slow)
	expect_true(all(abs(r$fast[fin] - r$slow[fin]) < 1e-4))   # < 0.1 mm

	r <- dist_both_v(rv, line, method = mth)
	expect_equal(is.na(r$fast), is.na(r$slow))
	fin <- !is.na(r$fast) & !is.na(r$slow)
	expect_true(all(abs(r$fast[fin] - r$slow[fin]) < 1e-4))
}

# unit = "km" branch for geo vector fast path
r <- dist_both_v(rv, poly, method = "geo", unit = "km")
expect_equal(is.na(r$fast), is.na(r$slow))
fin <- !is.na(r$fast) & !is.na(r$slow)
expect_true(all(abs(r$fast[fin] - r$slow[fin]) < 1e-7))
