
# planar hex grid: exact equal-area tessellation
e <- ext(0, 100, 0, 100)
h <- tessellate(e, size = 10)
expect_true(nrow(h) > 0)
a <- expanse(h, transform = FALSE)
# all planar cells have identical Cartesian area
expect_true((max(a) - min(a)) / mean(a) < 1e-9)


# flat-top variant
h2 <- tessellate(e, size = 12, flat_top = TRUE)
expect_true(nrow(h2) > 0)
a2 <- expanse(h2, transform = FALSE)
expect_true((max(a2) - min(a2)) / mean(a2) < 1e-9)


# global geo=TRUE: equal-area on the sphere
g <- tessellate(ext(-180, 180, -85, 85), size = 500000, geo = TRUE)
expect_true(nrow(g) == 2320)
ag <- expanse(g)
# spherical equal-area cells; spread should be small (< 2%)
expect_true(diff(range(ag)) / mean(ag) < 0.02)

# adjacent regional grids share the global anchor: cells whose
# centroids fall in the overlap region must coincide exactly
g_eu <- tessellate(ext(-10, 40, 30, 70), size = 200000, geo = TRUE)
g_af <- tessellate(ext( 10, 60, -5, 35), size = 200000, geo = TRUE)
cn_eu <- crds(centroids(g_eu, inside = FALSE))
cn_af <- crds(centroids(g_af, inside = FALSE))
in_overlap <- function(cn) {
	cn[, 1] >= 10.5 & cn[, 1] <= 39.5 & cn[, 2] >= 30.5 & cn[, 2] <= 34.5
}
sub_eu <- cn_eu[in_overlap(cn_eu), , drop = FALSE]
sub_af <- cn_af[in_overlap(cn_af), , drop = FALSE]
expect_equal(nrow(sub_eu), nrow(sub_af))
oe <- order(sub_eu[, 1], sub_eu[, 2])
oa <- order(sub_af[, 1], sub_af[, 2])
expect_equal(sub_eu[oe, , drop = FALSE], sub_af[oa, , drop = FALSE])


# antimeridian: cells that straddle +/-180 are returned as multipart
# geometries spanning lon -180 .. +180, with no duplicates and exactly
# Nx unique cells per row
h180 <- tessellate(ext(-180, 180, -80, 80), size = 500000, geo = TRUE)
cn <- crds(centroids(h180, inside = FALSE))
# count cells per latitude row (rounded centroid lat)
lat_round <- round(cn[, 2], 4)
counts <- as.integer(table(lat_round))
expect_true(all(counts == counts[1]))

# input validation
expect_error(tessellate(e, size = 0))
expect_error(tessellate(e, size = -1))
expect_error(tessellate(e, size = c(1, 2)))
expect_error(tessellate(n = 0, type="polyhedron"))
expect_error(tessellate(size = 0, type="polyh"))

# truncated dodecahedron-like base case (Goldberg G(2, 0))
g <- tessellate(n = 2, type="polyhedron")
expect_equal(nrow(g), 42)                          # 10 n^2 + 2
expect_equal(sum(g$type == "pentagon"), 12L)


# n = 10 -> 1002 cells, 12 pentagons + 990 hexagons
g10 <- tessellate(n=10, type="polyhedron")
expect_equal(nrow(g10), 1002)
expect_equal(sum(g10$type == "pentagon"), 12L)
expect_equal(sum(g10$type == "hexagon"),  990L)
expect_true(all(is.valid(g10)))


# cell areas should sum, at every resolution, to exactly the same total
ref_area <- sum(expanse(tessellate(n = 2, type="polyhedron")))
for (n in c(2L, 5L, 10L, 20L)) {
	g <- tessellate(n = n, type="polyhedron")
	expect_equal(nrow(g), 10L * n * n + 2)
	expect_true(all(is.valid(g)))
	a <- sum(expanse(g))
	expect_true(abs(a - ref_area) / ref_area < 1e-6)
}


# two of the 12 pentagons must sit on the geographic poles (their
# extent must reach lat = +/- 90).
g <- tessellate(n = 5, type="polyhedron")
pent <- g[g$type == "pentagon", ]
ymax <- vapply(seq_len(nrow(pent)), function(i) ymax(ext(pent[i, ])), numeric(1))
ymin <- vapply(seq_len(nrow(pent)), function(i) ymin(ext(pent[i, ])), numeric(1))
expect_true(any(ymax >=  90 - 1e-9))
expect_true(any(ymin <= -90 + 1e-9))


# size argument: derives n inversely from cell size
g_size <- tessellate(size = 2e6, type="polyhedron")
expect_true(nrow(g_size) > 0)
expect_true(all(g_size$type %in% c("pentagon", "hexagon")))
# smaller size -> more cells
g_small <- tessellate(size = 5e5, type="polyhedron")
expect_true(nrow(g_small) > nrow(g_size))


# extent filter selects a subset, including the polar extent edge cases
g_eq <- tessellate(ext(-30, 30, -10, 10), n = 10, type="polyhedron", geo=TRUE)
expect_true(nrow(g_eq) > 0)
expect_true(nrow(g_eq) < nrow(g10))


