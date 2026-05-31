# Author: Robert J. Hijmans
# Date :  2026
# Version 1.0
# License GPL v3


# Build a SpatNetwork from a SpatVector of lines.
setMethod("as.network", signature(x="SpatVector"),
	function(x, snap=0, merge=TRUE) {
		gt <- geomtype(x)
		if (gt != "lines") {
			error("as.network", "x must be a SpatVector of lines")
		}
		if (length(snap) != 1 || !is.numeric(snap) || is.na(snap) || snap < 0) {
			error("as.network", "snap must be a single non-negative number")
		}
		n <- methods::new("SpatNetwork")
		n@pntr <- x@pntr$as_network(snap, isTRUE(merge))
		messages(n, "as.network")
	}
)


setMethod("show", signature(object="SpatNetwork"),
	function(object) {
		cat(object@pntr$show())
	}
)


setMethod("nnodes", signature(x="SpatNetwork"),
	function(x, ...) {
		x@pntr$nnodes()
	}
)


setMethod("nedges", signature(x="SpatNetwork"),
	function(x, ...) {
		x@pntr$nedges()
	}
)


setMethod("nodes", signature(x="SpatNetwork"),
	function(x, ...) {
		v <- methods::new("SpatVector")
		v@pntr <- x@pntr$as_nodes()
		messages(v, "nodes")
	}
)


setMethod("edges", signature(x="SpatNetwork"),
	function(x, ...) {
		v <- methods::new("SpatVector")
		v@pntr <- x@pntr$as_edges()
		messages(v, "edges")
	}
)


setMethod("crs", signature(x="SpatNetwork"),
	function(x, proj=FALSE, describe=FALSE, parse=FALSE) {
		if (proj) {
			x@pntr$getSRS("proj4")
		} else {
			x@pntr$getSRS("wkt")
		}
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
		e <- edges(x)
		plot(e, col=col, ...)
		if (nodes) {
			points(nodes(x), col=node.col, cex=cex, pch=node.pch)
		}
		invisible(x)
	}
)
