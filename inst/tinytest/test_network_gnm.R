# GNM (GDAL Geographic Network Model) round-trip tests for SpatNetwork.
# Exercises writeNetwork() (write side) and netw(<filename>) (read side).

if (terra::gdal() >= "3.1.0" && "GNMFile" %in% terra::gdal(drivers = TRUE)$name) {

	# --- a small undirected, weighted test network --------------------
	a <- vect("LINESTRING(0 0, 10 10)", crs = "EPSG:32633")
	b <- vect("LINESTRING(0 10, 10 0)", crs = "EPSG:32633")
	roads <- rbind(a, b)
	roads$name <- c("A", "B")
	net <- netw(roads)

	# Round-trip: write, read, and compare.
	dst <- file.path(tempdir(), paste0("gnm_test_", as.integer(runif(1)*1e9)))
	on.exit(unlink(dst, recursive = TRUE, force = TRUE), add = TRUE)
	expect_silent(ok <- writeNetwork(net, dst))
	expect_true(isTRUE(ok))
	expect_true(dir.exists(dst))

	# GNMFile uses ESRI Shapefile under the hood by default; the system
	# layers (_gnm_meta, _gnm_graph, _gnm_features, _gnm_srs.prj) and
	# the class layers (nodes, edges) should all be present.
	files <- list.files(dst)
	expect_true("_gnm_graph.dbf" %in% files)
	expect_true("_gnm_meta.dbf"  %in% files)
	expect_true("_gnm_srs.prj"   %in% files)
	expect_true("nodes.shp"      %in% files)
	expect_true("edges.shp"      %in% files)

	# `netw(<filename>)` is the (only) read entry point.
	net2 <- netw(dst)
	expect_inherits(net2, "SpatNetwork")
	expect_equal(net_nnodes(net2), net_nnodes(net))
	expect_equal(net_nedges(net2), net_nedges(net))
	expect_equal(net_directed(net2), net_directed(net))
	expect_equal(sort(net_weights(net2)), sort(net_weights(net)))

	# CRS preserved (verbatim WKT may not match, but the EPSG-recognised
	# representation should round-trip).
	expect_true(nzchar(crs(net2)))
	expect_equal(crs(net2, describe = TRUE)$code,
	             crs(net,  describe = TRUE)$code)

	# Shortest paths should agree.
	sp1 <- shortestPath(net,  from = 1, to = 4)$distance
	sp2 <- shortestPath(net2, from = 1, to = 4)$distance
	expect_equal(sp1, sp2)

	# --- directed network --------------------------------------------
	ndir <- netw(roads, directed = TRUE)
	dst2 <- file.path(tempdir(), paste0("gnm_dir_", as.integer(runif(1)*1e9)))
	on.exit(unlink(dst2, recursive = TRUE, force = TRUE), add = TRUE)
	expect_silent(writeNetwork(ndir, dst2))
	ndir2 <- netw(dst2)
	expect_true(net_directed(ndir2))

	# --- error paths --------------------------------------------------

	# writing twice without overwrite must error
	expect_error(writeNetwork(net, dst))
	# but with overwrite=TRUE works
	expect_true(isTRUE(writeNetwork(net, dst, overwrite = TRUE)))

	# reading a missing path errors
	expect_error(netw(file.path(tempdir(), "definitely_not_a_gnm")))

	# bogus filetype is rejected
	expect_error(writeNetwork(net, dst, filetype = "Bogus", overwrite = TRUE))
}
