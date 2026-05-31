# SpatNetwork: build a planar topological network from lines.

# Two crossing roads -> 5 nodes, 4 edges
a <- vect("LINESTRING(0 0, 10 10)", crs = "local")
b <- vect("LINESTRING(0 10, 10 0)", crs = "local")
roads <- rbind(a, b)
roads$name <- c("A", "B")

net <- as.network(roads)

expect_inherits(net, "SpatNetwork")
expect_equal(nnodes(net), 5)
expect_equal(nedges(net), 4)

nds <- nodes(net)
eds <- edges(net)
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
expect_true(all(eds$from_node >= 1 & eds$from_node <= nnodes(net)))
expect_true(all(eds$to_node   >= 1 & eds$to_node   <= nnodes(net)))

# CRS is preserved end-to-end.
expect_equal(crs(nds), crs(roads))
expect_equal(crs(eds), crs(roads))


# Single line, no crossings -> 2 nodes, 1 edge, no warnings.
single <- vect("LINESTRING(0 0, 1 1, 2 0)", crs = "local")
net1 <- as.network(single)
expect_equal(nnodes(net1), 2)
expect_equal(nedges(net1), 1)
# Interior vertex (1,1) is preserved on the edge geometry.
e1 <- edges(net1)
expect_equal(nrow(crds(e1)), 3)


# T-junction: line C joins line A halfway along.
a2 <- vect("LINESTRING(0 0, 10 0)", crs = "local")
c2 <- vect("LINESTRING(5 0, 5 10)", crs = "local")
net2 <- as.network(rbind(a2, c2))
expect_equal(nnodes(net2), 4)  # endpoints (0,0),(10,0),(5,10) + junction (5,0)
expect_equal(nedges(net2), 3)  # A is split in two; C stays whole
junction_deg <- max(nodes(net2)$degree)
expect_equal(junction_deg, 3L)


# Wrong geometry type -> error.
pts <- vect(cbind(c(0, 1), c(0, 1)), crs = "local")
expect_error(as.network(pts))

# Negative snap -> error.
expect_error(as.network(roads, snap = -1))


# Degree-2 nodes from segmented inputs are removed by default.
# A single straight road digitized as three abutting segments should
# collapse to one edge with two endpoints under merge=TRUE.
seg1 <- vect("LINESTRING(0 0, 5 0)",  crs = "local")
seg2 <- vect("LINESTRING(5 0, 10 0)", crs = "local")
seg3 <- vect("LINESTRING(10 0, 15 0)", crs = "local")
chain <- rbind(seg1, seg2, seg3)

netM <- as.network(chain)              # merge = TRUE (default)
expect_equal(nnodes(netM), 2)
expect_equal(nedges(netM), 1)
# All non-extreme nodes should be gone, so degrees are exactly 1.
expect_equal(sort(nodes(netM)$degree), c(1L, 1L))

netU <- as.network(chain, merge = FALSE)
expect_equal(nnodes(netU), 4)
expect_equal(nedges(netU), 3)
# The two interior nodes are degree-2 join points.
expect_equal(sum(nodes(netU)$degree == 2L), 2L)

# A real junction is still preserved under merge=TRUE.
trunk <- vect("LINESTRING(0 0, 5 0, 10 0)", crs = "local")
spur  <- vect("LINESTRING(5 0, 5 5)",       crs = "local")
netJ <- as.network(rbind(trunk, spur))
expect_equal(nnodes(netJ), 4)  # (0,0) (5,0) (10,0) (5,5)
expect_equal(nedges(netJ), 3)
expect_equal(max(nodes(netJ)$degree), 3L)
