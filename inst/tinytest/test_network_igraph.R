# Coercion tests for SpatNetwork <-> igraph.

if (!requireNamespace("igraph", quietly = TRUE)) {
	exit_file("igraph not available")
}

a <- vect("LINESTRING(0 0, 10 10)", crs = "local")
b <- vect("LINESTRING(0 10, 10 0)", crs = "local")
roads <- rbind(a, b)
roads$name <- c("A", "B")
net <- netw(roads)


# --- SpatNetwork -> igraph -------------------------------------------------
g <- as(net, "igraph")
expect_inherits(g, "igraph")
expect_equal(igraph::vcount(g), net_nnodes(net))
expect_equal(igraph::ecount(g), net_nedges(net))
expect_false(igraph::is_directed(g))

# Vertex attributes carry the coordinates.
expect_true(all(c("x", "y") %in% igraph::vertex_attr_names(g)))
xy <- crds(net_nodes(net))
expect_equivalent(igraph::vertex_attr(g, "x"), xy[, 1])
expect_equivalent(igraph::vertex_attr(g, "y"), xy[, 2])

# Edge weights land on the conventional `weight` attribute.
expect_true("weight" %in% igraph::edge_attr_names(g))
expect_equivalent(igraph::edge_attr(g, "weight"), net_weights(net))

# `length` rides along.
expect_true("length" %in% igraph::edge_attr_names(g))
expect_equivalent(igraph::edge_attr(g, "length"), net_edges(net)$length)

# CRS round-trips through the `crs` graph attribute.
expect_true("crs" %in% igraph::graph_attr_names(g))


# Forward results agree with shortestPath: a corner-to-corner shortest path
# computed by Dijkstra in terra has the same total weight as igraph's.
sp_terra  <- shortestPath(net, from = 1, to = 4)$distance
sp_igraph <- as.numeric(igraph::distances(g, v = 1, to = 4))
expect_equal(sp_terra, sp_igraph, tolerance = 1e-9)


# Directed networks transfer their orientation.
ndir <- netw(roads, directed = TRUE)
gdir <- as(ndir, "igraph")
expect_true(igraph::is_directed(gdir))


# --- igraph -> SpatNetwork -------------------------------------------------
net2 <- as(g, "SpatNetwork")
expect_inherits(net2, "SpatNetwork")
expect_equal(net_nnodes(net2), net_nnodes(net))
expect_equal(net_nedges(net2), net_nedges(net))
expect_equal(net_directed(net2), net_directed(net))

# `netw(igraph)` is an equivalent entry point.
net2b <- netw(g)
expect_equal(net_nnodes(net2b), net_nnodes(net))
expect_equal(net_nedges(net2b), net_nedges(net))

# Node coordinates survive intact.
expect_equivalent(crds(net_nodes(net2)), crds(net_nodes(net)))

# Weights survive too (within tolerance: same numeric values, possibly
# reordered if igraph reorders internally -- which it doesn't here).
expect_equivalent(sort(net_weights(net2)), sort(net_weights(net)), tolerance = 1e-9)

# CRS round-trips through the graph attribute.
expect_equal(crs(net2), crs(net))


# Building from a fresh igraph (no terra origin): need x / y vertex attrs.
g2 <- igraph::make_graph(c(1,2, 2,3, 3,1), directed = FALSE)
g2 <- igraph::set_vertex_attr(g2, "x", value = c(0, 1, 0))
g2 <- igraph::set_vertex_attr(g2, "y", value = c(0, 0, 1))
g2 <- igraph::set_edge_attr(g2,   "weight", value = c(1, 2, 3))
nfresh <- netw(g2, crs = "local")
expect_equal(net_nnodes(nfresh), 3)
expect_equal(net_nedges(nfresh), 3)
expect_equivalent(sort(net_weights(nfresh)), c(1, 2, 3))

# Missing coordinate attributes -> error.
g3 <- igraph::make_graph(c(1,2, 2,3), directed = FALSE)
expect_error(netw(g3))
expect_error(as(g3, "SpatNetwork"))


# Lossless idempotency on already-a-SpatNetwork.
expect_identical(netw(net), net)


# Shortest paths agree on a converted graph.
sp_pre  <- shortestPath(net,  from = 1, to = 4)$distance
sp_post <- shortestPath(net2, from = 1, to = 4)$distance
expect_equal(sp_pre, sp_post, tolerance = 1e-9)
