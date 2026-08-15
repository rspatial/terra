## Z coordinate storage and round-trips (#824)

v2 <- vect(cbind(1:3, 10:12))
expect_false(has.z(v2))

v3 <- vect()
v3@pntr$setPointsXYZ(as.double(1:3), as.double(10:12), as.double(c(100, 200, 300)))
expect_true(has.z(v3))
expect_equal(nrow(v3), 3L)
xy <- crds(v3)
expect_equivalent(xy[,1], 1:3)
expect_equivalent(xy[,2], 10:12)
expect_equivalent(xy[,3], c(100, 200, 300))

## empty / wrong-length z falls back to 2D
v4 <- vect()
v4@pntr$setPointsXYZ(as.double(1:2), as.double(3:4), as.double(numeric(0)))
expect_false(has.z(v4))

## geom() <-> vect() round-trip
g <- geom(v3)
expect_true("z" %in% colnames(g))
expect_equivalent(g[,"z"], c(100, 200, 300))
v3b <- vect(g, type="points")
expect_true(has.z(v3b))
expect_equivalent(crds(v3b)[,"z"], c(100, 200, 300))

## writeVector / vect(file) round-trip (GPKG)
f <- tempfile(fileext=".gpkg")
writeVector(v3, f, overwrite=TRUE)
v3c <- vect(f)
expect_true(has.z(v3c))
expect_equivalent(crds(v3c)[,"z"], c(100, 200, 300))

## sf coercion round-trip (if sf is available)
if (requireNamespace("sf", quietly=TRUE)) {
	s <- sf::st_as_sf(data.frame(x=1:3, y=10:12, z=c(100,200,300)),
		coords=c("x","y","z"), crs="EPSG:4326")
	vs <- vect(s)
	expect_true(has.z(vs))
	expect_equivalent(crds(vs)[,"z"], c(100, 200, 300))

	# SpatVector -> sf -> SpatVector
	s2 <- sf::st_as_sf(v3)
	vs2 <- vect(s2)
	expect_true(has.z(vs2))
	expect_equivalent(crds(vs2)[,"z"], c(100, 200, 300))
}
