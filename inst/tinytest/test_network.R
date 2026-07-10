# SpatNetwork: build a planar topological network from lines.

# Two crossing roads -> 5 nodes, 4 edges
a <- vect("LINESTRING(0 0, 10 10)", crs = "local")
b <- vect("LINESTRING(0 10, 10 0)", crs = "local")
roads <- rbind(a, b)
roads$name <- c("A", "B")

net <- netw(roads)

expect_inherits(net, "SpatNetwork")
expect_equal(net_nnodes(net), 5)
expect_equal(net_nedges(net), 4)

nds <- net_nodes(net)
eds <- net_edges(net)
expect_inherits(nds, "SpatVector")
expect_inherits(eds, "SpatVector")
expect_equal(geomtype(nds), "points")
expect_equal(geomtype(eds), "lines")

# Crossing node has degree 4; the four endpoints have degree 1.
deg <- nds$degree
expect_equal(sort(deg), c(1L, 1L, 1L, 1L, 4L))

# Edge lengths are all sqrt(50) (half-diagonal of the 10x10 box).
expect_equivalent(eds$length, rep(sqrt(50), 4), tolerance = 1e-9)

# Source attribution: 1-based index pointing back into the input,
# and the per-source attribute (`name`) is forwarded.
expect_true(all(eds$source_id %in% c(1L, 2L)))
expect_true(all(eds$name %in% c("A", "B")))
expect_equal(eds$name, roads$name[eds$source_id])

# Node ordering: from_node / to_node are valid 1-based indices.
expect_true(all(eds$from_node >= 1 & eds$from_node <= net_nnodes(net)))
expect_true(all(eds$to_node   >= 1 & eds$to_node   <= net_nnodes(net)))

# CRS is preserved end-to-end.
expect_equal(crs(nds), crs(roads))
expect_equal(crs(eds), crs(roads))

# Coercion via setAs gives the same network.
net_via_as <- as(roads, "SpatNetwork")
expect_inherits(net_via_as, "SpatNetwork")
expect_equal(net_nnodes(net_via_as), net_nnodes(net))
expect_equal(net_nedges(net_via_as), net_nedges(net))


# Single line, no crossings -> 2 nodes, 1 edge, no warnings.
single <- vect("LINESTRING(0 0, 1 1, 2 0)", crs = "local")
net1 <- netw(single)
expect_equal(net_nnodes(net1), 2)
expect_equal(net_nedges(net1), 1)
# Interior vertex (1,1) is preserved on the edge geometry.
e1 <- net_edges(net1)
expect_equal(nrow(crds(e1)), 3)


# T-junction: line C joins line A halfway along.
a2 <- vect("LINESTRING(0 0, 10 0)", crs = "local")
c2 <- vect("LINESTRING(5 0, 5 10)", crs = "local")
net2 <- netw(rbind(a2, c2))
expect_equal(net_nnodes(net2), 4)  # endpoints (0,0),(10,0),(5,10) + junction (5,0)
expect_equal(net_nedges(net2), 3)  # A is split in two; C stays whole
junction_deg <- max(net_nodes(net2)$degree)
expect_equal(junction_deg, 3L)


# Wrong geometry type -> error.
pts <- vect(cbind(c(0, 1), c(0, 1)), crs = "local")
expect_error(netw(pts))

# Negative snap -> error.
expect_error(netw(roads, snap = -1))


# Degree-2 nodes from segmented inputs are removed by default.
# A single straight road digitized as three abutting segments should
# collapse to one edge with two endpoints under merge=TRUE.
seg1 <- vect("LINESTRING(0 0, 5 0)",  crs = "local")
seg2 <- vect("LINESTRING(5 0, 10 0)", crs = "local")
seg3 <- vect("LINESTRING(10 0, 15 0)", crs = "local")
chain <- rbind(seg1, seg2, seg3)

netM <- netw(chain)              # merge = TRUE (default)
expect_equal(net_nnodes(netM), 2)
expect_equal(net_nedges(netM), 1)
# All non-extreme nodes should be gone, so degrees are exactly 1.
expect_equal(sort(net_nodes(netM)$degree), c(1L, 1L))

netU <- netw(chain, merge = FALSE)
expect_equal(net_nnodes(netU), 4)
expect_equal(net_nedges(netU), 3)
# The two interior nodes are degree-2 join points.
expect_equal(sum(net_nodes(netU)$degree == 2L), 2L)

# A real junction is still preserved under merge=TRUE.
trunk <- vect("LINESTRING(0 0, 5 0, 10 0)", crs = "local")
spur  <- vect("LINESTRING(5 0, 5 5)",       crs = "local")
netJ <- netw(rbind(trunk, spur))
expect_equal(net_nnodes(netJ), 4)  # (0,0) (5,0) (10,0) (5,5)
expect_equal(net_nedges(netJ), 3)
expect_equal(max(net_nodes(netJ)$degree), 3L)


# --- directed / weighted ---
# Default: undirected, weighted by length.
expect_false(net_directed(net))
w <- net_weights(net)
expect_true(!is.null(w))
expect_equal(length(w), net_nedges(net))
expect_equivalent(w, rep(sqrt(50), 4), tolerance = 1e-9)

# Explicitly unweighted -> net_weights() == NULL.
netU <- netw(roads, weights = FALSE)
expect_null(net_weights(netU))

# Custom numeric weights via constructor.
netCW <- netw(roads, weights = c(1, 2, 3, 4))
expect_equal(net_weights(netCW), c(1, 2, 3, 4))

# Custom weights via assignment.
netW <- netw(roads)
net_weights(netW) <- 10 * net_weights(netW)
expect_equivalent(net_weights(netW), 10 * sqrt(50) * rep(1, 4), tolerance = 1e-9)
# Wrong length is rejected.
expect_error(net_weights(netW) <- 1:3)
# Setting NULL clears weights.
net_weights(netW) <- NULL
expect_null(net_weights(netW))

# Directed networks expose in_degree / out_degree on net_nodes().
netD <- netw(roads, directed = TRUE)
expect_true(net_directed(netD))
nd <- net_nodes(netD)
expect_true(all(c("degree", "in_degree", "out_degree") %in% names(nd)))
# Sum of in_degree across nodes must equal nedges.
expect_equal(sum(nd$in_degree),  net_nedges(netD))
expect_equal(sum(nd$out_degree), net_nedges(netD))

# net_edges() exposes a "weight" column when weighted, not when unweighted.
expect_true("weight" %in% names(net_edges(net)))
expect_false("weight" %in% names(net_edges(netU)))

# Lon/lat networks: length is in meters via geodesic distance.
gll <- vect("LINESTRING(0 0, 1 0)", crs = "EPSG:4326")
netll <- netw(gll)
# 1 degree of longitude at the equator on WGS84 ~= 111319.49 m.
expect_equal(net_weights(netll), 111319.49, tolerance = 1)
expect_equal(net_edges(netll)$length, 111319.49, tolerance = 1)


# --- weights from a source column ---
roads2 <- rbind(a, b)
roads2$cost <- c(7, 100)
netCol <- netw(roads2, weights = "cost")
# Each input feature should have produced 2 edges, so the weight
# vector contains exactly two 7s and two 100s.
expect_equal(sort(net_weights(netCol)), c(7, 7, 100, 100))
# Unknown column -> error.
expect_error(netw(roads2, weights = "no_such_column"))
# Non-numeric column -> error.
roads2$nm <- c("a", "b")
expect_error(netw(roads2, weights = "nm"))


# --- shortestPath ---
# Two crossing roads, same length, undirected, weighted by length.
sp <- shortestPath(net, from = 1, to = 4)
expect_inherits(sp, "SpatVector")
expect_equal(nrow(sp), 1L)
expect_equal(sp$from, 1L)
expect_equal(sp$to,   4L)
# Path goes through the centre node and is exactly two edges.
expect_equal(sp$distance, 2 * sqrt(50), tolerance = 1e-9)
# Geometry is non-empty and traverses the centre point.
gxy <- geom(sp)[, c("x","y")]
expect_true(any(abs(gxy[,1] - 5) < 1e-9 & abs(gxy[,2] - 5) < 1e-9))

# from == to: distance 0, empty geometry.
sp0 <- shortestPath(net, from = 2, to = 2)
expect_equal(sp0$distance, 0)
expect_equal(nrow(crds(sp0)), 0L)

# Multiple destinations from one source -> recycled.
spM <- shortestPath(net, from = 1, to = 1:5)
expect_equal(nrow(spM), 5L)
expect_equal(spM$from, rep(1L, 5))
expect_equal(spM$to,   1:5)
expect_equal(spM$distance[1], 0)
# Distance from corner to opposite corner == 2*sqrt(50).
expect_equal(spM$distance[4], 2 * sqrt(50), tolerance = 1e-9)

# Custom weights via a column propagate into shortestPath.
spC <- shortestPath(netCol, from = 1, to = 4)
# The shortest cost path uses the cheap (cost=7) road -> has to reroute.
# At minimum the distance must be <= the all-expensive route.
expect_true(spC$distance <= 2 * 100 + 1e-9)

# Disconnected network -> NA distance and empty geom.
disc1 <- vect("LINESTRING(0 0, 1 1)", crs = "local")
disc2 <- vect("LINESTRING(10 10, 11 11)", crs = "local")
ndisc <- netw(rbind(disc1, disc2))
spD <- shortestPath(ndisc, from = 1, to = net_nnodes(ndisc))
expect_true(is.na(spD$distance))
expect_equal(nrow(crds(spD)), 0L)

# Directed networks respect orientation.
# Build a Y: trunk 1 -> 2 -> 3 with a spur 2 -> 4.
trunk2 <- vect("LINESTRING(0 0, 5 0, 10 0)", crs = "local")
spur2  <- vect("LINESTRING(5 0, 5 5)",       crs = "local")
ndir   <- netw(rbind(trunk2, spur2), directed = TRUE)
# Forward should work.
spF <- shortestPath(ndir, from = 1, to = 3)
expect_false(is.na(spF$distance))
# Reverse should be unreachable in a directed graph.
spR <- shortestPath(ndir, from = 3, to = 1)
expect_true(is.na(spR$distance))
# Same on the undirected version: both directions reachable.
nund <- netw(rbind(trunk2, spur2), directed = FALSE)
expect_false(is.na(shortestPath(nund, from = 3, to = 1)$distance))

# Invalid node ids -> error.
expect_error(shortestPath(net, from = 0,    to = 1))
expect_error(shortestPath(net, from = 9999, to = 1))


# --- empty network from netw() ---
empty <- netw()
expect_inherits(empty, "SpatNetwork")
expect_equal(net_nnodes(empty), 0)
expect_equal(net_nedges(empty), 0)


# --- netw() is idempotent on SpatNetwork ---
expect_identical(netw(net), net)
