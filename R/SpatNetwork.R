# Author: Robert J. Hijmans
# Date :  2026
# Version 1.0
# License GPL v3


# ---------------------------------------------------------------------------
# Internal builders -- not exported. All user-facing entry points
# (`netw()`, `setAs(..., "SpatNetwork")`) ultimately call these.
# ---------------------------------------------------------------------------

.netw_from_SpatVector <- function(x, snap=0, merge=TRUE, directed=FALSE, weights=TRUE) {
	gt <- geomtype(x)
	if (gt != "lines") {
		error("netw", "x must be a SpatVector of lines")
	}
	if (length(snap) != 1 || !is.numeric(snap) || is.na(snap) || snap < 0) {
		error("netw", "snap must be a single non-negative number")
	}

	# `weights` argument: TRUE/FALSE, a numeric vector (length nedges),
	# or a single character matching a numeric column of `x` (one
	# weight per source feature, looked up by source_id after noding).
	w_user_vec <- NULL
	w_colname  <- NULL
	w_flag     <- FALSE
	if (is.logical(weights)) {
		if (length(weights) != 1 || is.na(weights)) {
			error("netw", "'weights' must be TRUE/FALSE, a numeric vector, or a column name")
		}
		w_flag <- isTRUE(weights)
	} else if (is.numeric(weights)) {
		w_user_vec <- as.numeric(weights)
	} else if (is.character(weights) && length(weights) == 1 && !is.na(weights) && nzchar(weights)) {
		if (!(weights %in% names(x))) {
			error("netw", paste0("'", weights, "' is not a column of x"))
		}
		col <- x[[weights, drop=TRUE]]
		if (!is.numeric(col)) {
			error("netw", paste0("column '", weights, "' must be numeric"))
		}
		w_colname <- weights
		w_flag    <- TRUE
	} else {
		error("netw", "'weights' must be TRUE/FALSE, a numeric vector, or a column name")
	}

	n <- methods::new("SpatNetwork")
	n@pntr <- x@pntr$as_network(snap, isTRUE(merge), isTRUE(directed), w_flag)
	n <- messages(n, "netw")

	if (!is.null(w_user_vec)) {
		if (length(w_user_vec) != net_nedges(n)) {
			error("netw", "length of 'weights' must equal the number of edges")
		}
		n@pntr$setWeights(w_user_vec)
		n <- messages(n, "netw")
	} else if (!is.null(w_colname)) {
		# Look up per-edge weight via source_id (1-based; -1 = unattributed).
		eds <- net_edges(n)
		src <- eds$source_id
		col <- x[[w_colname, drop=TRUE]]
		w   <- rep(NA_real_, length(src))
		ok  <- src >= 1L & !is.na(src)
		w[ok] <- col[src[ok]]
		bad <- which(is.na(w))
		if (length(bad) > 0) {
			len <- eds$length
			w[bad] <- len[bad]
			warn("netw", paste("falling back to edge length for",
				length(bad), "edge(s) with no source attribution"))
		}
		n@pntr$setWeights(as.numeric(w))
		n <- messages(n, "netw")
	}
	n
}


.netw_from_igraph <- function(x, x_attr = "x", y_attr = "y",
                              weight_attr = "weight", crs = NULL) {
	if (!requireNamespace("igraph", quietly = TRUE)) {
		error("netw", "the 'igraph' package is required; install it first")
	}
	if (!(x_attr %in% igraph::vertex_attr_names(x)) ||
	    !(y_attr %in% igraph::vertex_attr_names(x))) {
		error("netw", paste0("the igraph object must have vertex attributes '",
			x_attr, "' and '", y_attr, "' giving node coordinates"))
	}
	nx <- as.numeric(igraph::vertex_attr(x, x_attr))
	ny <- as.numeric(igraph::vertex_attr(x, y_attr))
	if (anyNA(nx) || anyNA(ny)) {
		error("netw", "igraph node coordinates contain NA")
	}

	el <- igraph::as_edgelist(x, names = FALSE)   # 1-based vertex ids
	if (NROW(el) == 0) {
		efrom <- integer(0)
		eto   <- integer(0)
	} else {
		efrom <- as.integer(el[, 1])
		eto   <- as.integer(el[, 2])
	}

	w <- numeric(0)
	if (weight_attr %in% igraph::edge_attr_names(x)) {
		w <- as.numeric(igraph::edge_attr(x, weight_attr))
	}
	directed <- igraph::is_directed(x)

	n <- methods::new("SpatNetwork")
	n@pntr <- SpatNetwork$new()
	ok <- n@pntr$buildFromComponents(nx, ny,
	                                 as.integer(efrom - 1L),
	                                 as.integer(eto - 1L),
	                                 w, directed)
	n <- messages(n, "netw")
	if (!isTRUE(ok)) return(n)

	# CRS: prefer explicit `crs=`, then graph_attr "crs", then leave empty.
	# `.txtCRS` normalises shortcuts like "local" / EPSG codes.
	if (!is.null(crs)) {
		n@pntr$setSRS(.txtCRS(crs, warn = FALSE))
	} else if ("crs" %in% igraph::graph_attr_names(x)) {
		v <- as.character(igraph::graph_attr(x, "crs"))
		if (nzchar(v)) n@pntr$setSRS(.txtCRS(v, warn = FALSE))
	}
	messages(n, "netw")
}


.netw_to_igraph <- function(x) {
	if (!requireNamespace("igraph", quietly = TRUE)) {
		error("as,SpatNetwork,igraph", "the 'igraph' package is required; install it first")
	}
	eds <- net_edges(x)
	nds <- net_nodes(x)

	# Edge list as a 2-column integer matrix (1-based vertex ids).
	el <- cbind(eds$from_node, eds$to_node)
	g <- igraph::graph_from_edgelist(el, directed = net_directed(x))

	# Vertex attributes: coordinates plus the columns of net_nodes(x).
	xy <- crds(nds)
	g <- igraph::set_vertex_attr(g, "x", value = xy[, 1])
	g <- igraph::set_vertex_attr(g, "y", value = xy[, 2])
	for (nm in setdiff(names(nds), c("x", "y"))) {
		g <- igraph::set_vertex_attr(g, nm, value = nds[[nm, drop = TRUE]])
	}

	# Edge attributes: `weight` (when weighted) plus everything carried
	# on net_edges() except the structural columns.
	w <- net_weights(x)
	if (!is.null(w)) {
		g <- igraph::set_edge_attr(g, "weight", value = w)
	}
	for (nm in setdiff(names(eds), c("from_node", "to_node"))) {
		g <- igraph::set_edge_attr(g, nm, value = eds[[nm, drop = TRUE]])
	}

	# CRS travels along as a graph attribute so a round-trip preserves it.
	wkt <- crs(x)
	if (!is.na(wkt) && nzchar(wkt)) {
		g <- igraph::set_graph_attr(g, "crs", wkt)
	}
	g
}


# ---------------------------------------------------------------------------
# `netw()` -- the user-facing constructor, in the spirit of `vect()`.
#
# Dispatches by class:
#   netw()                      -> empty SpatNetwork
#   netw(SpatVector, ...)       -> build from line geometry (planar noding)
#   netw(igraph,     ...)       -> build from an igraph object
#   netw(<filename>, ...)       -> read a GDAL GNM dataset from disk
#   netw(SpatNetwork, ...)      -> identity (idempotent)
# ---------------------------------------------------------------------------

setMethod("netw", signature(x="missing"),
	function(x, ...) {
		n <- methods::new("SpatNetwork")
		n@pntr <- SpatNetwork$new()
		n
	}
)

setMethod("netw", signature(x="SpatVector"),
	function(x, snap=0, merge=TRUE, directed=FALSE, weights=TRUE, ...) {
		.netw_from_SpatVector(x, snap=snap, merge=merge,
		                      directed=directed, weights=weights)
	}
)

setMethod("netw", signature(x="igraph"),
	function(x, x_attr="x", y_attr="y", weight_attr="weight", crs=NULL, ...) {
		.netw_from_igraph(x, x_attr=x_attr, y_attr=y_attr,
		                  weight_attr=weight_attr, crs=crs)
	}
)

setMethod("netw", signature(x="character"),
	function(x, ...) {
		.gnm_read(x)
	}
)

setMethod("netw", signature(x="SpatNetwork"),
	function(x, ...) x
)


# ---------------------------------------------------------------------------
# Formal coercion via setAs (so `as(x, "SpatNetwork")` and
# `as(net, "igraph")` work, mirroring how other Spat* classes coerce).
# ---------------------------------------------------------------------------

setAs("SpatVector", "SpatNetwork",
	function(from) .netw_from_SpatVector(from)
)

setAs("igraph", "SpatNetwork",
	function(from) .netw_from_igraph(from)
)

setAs("SpatNetwork", "igraph",
	function(from) .netw_to_igraph(from)
)


# ---------------------------------------------------------------------------
# show / dimensions
# ---------------------------------------------------------------------------

setMethod("show", signature(object="SpatNetwork"),
	function(object) {
		cat(object@pntr$show())
	}
)

setMethod("net_nnodes", signature(x="SpatNetwork"),
	function(x, ...) {
		x@pntr$nnodes()
	}
)

setMethod("net_nedges", signature(x="SpatNetwork"),
	function(x, ...) {
		x@pntr$nedges()
	}
)


# ---------------------------------------------------------------------------
# nodes / edges (returned as SpatVectors of points / lines).
# Names use the `net_*` prefix to avoid collisions with igraph and the
# `network` package, both of which export plain `nodes`, `edges`, etc.
# ---------------------------------------------------------------------------

setMethod("net_nodes", signature(x="SpatNetwork"),
	function(x, ...) {
		v <- methods::new("SpatVector")
		v@pntr <- x@pntr$as_nodes()
		messages(v, "net_nodes")
	}
)

setMethod("net_edges", signature(x="SpatNetwork"),
	function(x, ...) {
		v <- methods::new("SpatVector")
		v@pntr <- x@pntr$as_edges()
		v <- messages(v, "net_edges")
		# C++ uses 0-based node indices, with -1 marking an unattributed
		# source feature. Translate to R-idiomatic 1-based ids and NA.
		v$from_node <- v$from_node + 1L
		v$to_node   <- v$to_node   + 1L
		sid <- v$source_id
		sid[sid < 0L] <- NA_integer_
		v$source_id <- sid + 1L
		v
	}
)


# ---------------------------------------------------------------------------
# CRS, extent, plotting
# ---------------------------------------------------------------------------

setMethod("crs", signature(x="SpatNetwork"),
	function(x, proj=FALSE, describe=FALSE, parse=FALSE) {
		wkt <- x@pntr$getSRS("wkt")
		if (describe) {
			# Reuse the SpatVector path so describe/parse behave the same
			# way as for other Spat* classes.
			v <- vect(cbind(0, 0), crs = wkt)
			return(.get_CRS(v, proj=proj, describe=describe, parse=parse))
		}
		if (proj) {
			r <- x@pntr$getSRS("proj4")
		} else {
			r <- wkt
		}
		if (parse) unlist(strsplit(r, "\n")) else r
	}
)

setMethod("ext", signature(x="SpatNetwork"),
	function(x, ...) {
		e <- methods::new("SpatExtent")
		e@pntr <- x@pntr$extent
		e
	}
)

setMethod("plot", signature(x="SpatNetwork", y="missing"),
	function(x, y, nodes=TRUE, col="black", node.col="red", node.pch=20, cex=1, ...) {
		e <- net_edges(x)
		plot(e, col=col, ...)
		if (nodes) {
			points(net_nodes(x), col=node.col, cex=cex, pch=node.pch)
		}
		invisible(x)
	}
)


# ---------------------------------------------------------------------------
# Directionality and weights. Renamed to net_* to avoid masking by
# igraph (which exports `is.directed`, `weights`, `nodes`, `edges`).
# ---------------------------------------------------------------------------

setMethod("net_directed", signature(x="SpatNetwork"),
	function(x, ...) {
		x@pntr$isDirected()
	}
)

setMethod("net_weights", signature(x="SpatNetwork"),
	function(x, ...) {
		if (!x@pntr$isWeighted()) return(NULL)
		x@pntr$getWeights()
	}
)

setReplaceMethod("net_weights", signature(x="SpatNetwork"),
	function(x, value) {
		if (is.null(value)) {
			x@pntr$clearWeights()
		} else {
			value <- as.numeric(value)
			if (length(value) != net_nedges(x)) {
				error("net_weights<-", "length of value must equal the number of edges")
			}
			x@pntr$setWeights(value)
		}
		messages(x, "net_weights<-")
	}
)


# ---------------------------------------------------------------------------
# Shortest paths between pairs of node ids.
#   from, to : integer node ids (1-based). Either may be length 1 to be
#              recycled against the other; otherwise lengths must match.
# Returns a SpatVector of lines, one feature per pair, with columns
# `from`, `to`, and `distance` (the total accumulated weight; NA when
# unreachable, 0 with empty geometry when from==to).
# ---------------------------------------------------------------------------

setMethod("shortestPath", signature(x="SpatNetwork"),
	function(x, from, to, ...) {
		if (missing(from) || missing(to)) {
			error("shortestPath", "both 'from' and 'to' must be provided")
		}
		from <- as.integer(from)
		to   <- as.integer(to)
		if (length(from) == 0 || length(to) == 0) {
			error("shortestPath", "'from' and 'to' must each have length >= 1")
		}
		if (anyNA(from) || anyNA(to)) {
			error("shortestPath", "'from' and 'to' must be integer node ids (no NAs)")
		}
		nn <- net_nnodes(x)
		if (any(from < 1L) || any(from > nn)) {
			error("shortestPath", "'from' contains an invalid node id")
		}
		if (any(to < 1L) || any(to > nn)) {
			error("shortestPath", "'to' contains an invalid node id")
		}

		v <- methods::new("SpatVector")
		# C++ side uses 0-based indices throughout.
		v@pntr <- x@pntr$shortest_paths(from - 1L, to - 1L)
		v <- messages(v, "shortestPath")
		v$from <- v$from + 1L
		v$to   <- v$to   + 1L
		v
	}
)


# ---------------------------------------------------------------------------
# Read / write the GDAL Geographic Network Model (GNM) on-disk format.
#
# `writeNetwork(net, filename, ...)` is the user-facing writer; for now it
# only knows about GNM (filetype "GNMFile" -- a directory of small
# OGR-format files -- or "GNMDatabase"), but the signature is structured
# so additional formats can be added later by dispatching on `filetype`.
#
# Reading is done through `netw(<filename>)`, which calls the internal
# `.gnm_read()`. We deliberately do not export a separate `gnm_read()`:
# users go through `netw()` for symmetry with `vect()` / `rast()`.
#
# A GNM dataset is a directory (driver "GNMFile") or a database
# (driver "GNMDatabase") containing two "class" layers (`nodes` and
# `edges`) plus the GNM system layers (`_gnm_meta`, `_gnm_graph`,
# `_gnm_features`, `_gnm_srs.prj`) that record topology, GFIDs and the
# network's spatial reference.
# ---------------------------------------------------------------------------

setMethod("writeNetwork", signature(x="SpatNetwork", filename="character"),
	function(x, filename, filetype="GNMFile", overwrite=FALSE, options=NULL, ...) {
		filename <- as.character(filename)[1]
		if (!nzchar(filename)) {
			error("writeNetwork", "'filename' must be a non-empty path")
		}
		filetype <- as.character(filetype)[1]
		if (!(filetype %in% c("GNMFile", "GNMDatabase"))) {
			error("writeNetwork",
			      paste0("unknown filetype '", filetype,
			             "'; supported: 'GNMFile', 'GNMDatabase'"))
		}
		options <- if (is.null(options)) character(0) else as.character(options)

		if (file.exists(filename)) {
			if (!isTRUE(overwrite)) {
				error("writeNetwork",
				      paste0("'", filename, "' already exists; ",
				             "set overwrite=TRUE to replace it"))
			}
			# GNM refuses to write into an existing path; remove first.
			unlink(filename, recursive = TRUE, force = TRUE)
		}

		ok <- x@pntr$write_gnm(filename, filetype, options)
		messages(x, "writeNetwork")
		invisible(isTRUE(ok))
	}
)


.gnm_read <- function(filename) {
	filename <- as.character(filename)[1]
	if (!nzchar(filename)) {
		error("netw", "'filename' must be a non-empty path")
	}
	if (!file.exists(filename)) {
		error("netw", paste0("'", filename, "' does not exist"))
	}
	n <- methods::new("SpatNetwork")
	n@pntr <- SpatNetwork$new()
	ok <- n@pntr$read_gnm(filename)
	n <- messages(n, "netw")
	if (!isTRUE(ok)) return(n)
	n
}
