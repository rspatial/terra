
#######################################################################
# Regressions reported via reverse-dependency checks (terra 1.9-22):
#   - SiMRiv / raptr : as(sf, "SpatVector") -> "non-conformable arrays"
#                      because the FULL-polygon 2x2 detection compared
#                      values before checking dimensions.
#   - mapmisc        : terra::project(SpatVector, ...) raised a hard
#                      error when all points became unprojectable
#                      (e.g. projecting the poles into a tilted-
#                      perspective view); should be a warning so that
#                      callers using suppressWarnings() can recover.
#######################################################################

# ------------------------------------------------------------------ #
# Fix 1: as(sf, "SpatVector") on an ordinary polygon must not error  #
#        with "non-conformable arrays"                                #
# ------------------------------------------------------------------ #
if (requireNamespace("sf", quietly = TRUE)) {

	# A polygon whose first ring has > 2 vertices: this is the case
	# that triggered the old "non-conformable arrays" error.
	pol <- sf::st_polygon(list(rbind(
		c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(0, 0)
	)))
	sfdf <- sf::st_sf(geometry = sf::st_sfc(pol, crs = 4326),
	                  id = 1L)

	v <- methods::as(sfdf, "SpatVector")
	expect_inherits(v, "SpatVector")
	expect_equal(nrow(v), 1L)
	expect_equal(geomtype(v), "polygons")

	# A second sanity check: also works inside rasterize()
	r <- rast(extent = ext(0, 1, 0, 1), nrows = 4, ncols = 4)
	rr <- rasterize(sfdf, r, field = "id")
	expect_inherits(rr, "SpatRaster")

	# ------------------------------------------------------------ #
	# sf encodes the FULL polygon (whole sphere) as a 2x2 matrix.   #
	# All three coercion paths -- sf, sfc, sfg -- must detect this  #
	# and return a global polygon (not error or compare values      #
	# blindly).                                                     #
	# ------------------------------------------------------------ #
	g_sfc <- sf::st_as_sfc("POLYGON FULL", crs = "EPSG:4326")

	# user-reported reproducer: vect(<sfc POLYGON FULL>)
	v_sfc <- vect(g_sfc)
	expect_inherits(v_sfc, "SpatVector")
	expect_equal(nrow(v_sfc), 1L)
	expect_equal(geomtype(v_sfc), "polygons")
	expect_equal(unname(as.vector(ext(v_sfc))), c(-180, 180, -90, 90))

	# wrapped in an sf data frame
	sf_full <- sf::st_sf(geometry = g_sfc, id = 1L)
	v_sf <- vect(sf_full)
	expect_inherits(v_sf, "SpatVector")
	expect_equal(nrow(v_sf), 1L)
	expect_equal(unname(as.vector(ext(v_sf))), c(-180, 180, -90, 90))

	# bare sfg (XY POLYGON)
	v_sfg <- vect(g_sfc[[1]])
	expect_inherits(v_sfg, "SpatVector")
	expect_equal(nrow(v_sfg), 1L)
	expect_equal(unname(as.vector(ext(v_sfg))), c(-180, 180, -90, 90))
}


# ------------------------------------------------------------------ #
# Fix 2: project() must not hard-error when all output coordinates   #
#        are unprojectable; emit a warning instead.                  #
# ------------------------------------------------------------------ #

# A tilted-perspective projection centered far from the poles
# is a setup that will produce NaN for the polar points.
tpers <- "+proj=tpers +h=2000000 +tilt=-15 +lat_0=15 +lon_0=130 +datum=WGS84"
poles <- vect(cbind(c(0, 0), c(90, -90)), crs = "EPSG:4326")

# Old behaviour: error stopping execution.
# New behaviour: warning, project() still returns a SpatVector.
projected <- suppressWarnings(project(poles, tpers))
expect_inherits(projected, "SpatVector")
expect_equal(nrow(projected), 2L)

# When wrapped in withCallingHandlers we should see a warning, not error.
saw_warning <- FALSE
res <- withCallingHandlers(
	project(poles, tpers),
	warning = function(w) {
		saw_warning <<- TRUE
		invokeRestart("muffleWarning")
	}
)
expect_inherits(res, "SpatVector")
expect_true(saw_warning)
