
setMethod("length", signature(x="SpatVectorCollection"),
	function(x) {
		x@pntr$size()
	}
)

setMethod("svc", signature(x="missing"),
	function(x) {
		v <- methods::new("SpatVectorCollection")
		v@pntr <- SpatVectorCollection$new()
		v
	}
)

setMethod("svc", signature(x="character"),
	function(x, layer="", query="", dialect="", extent=NULL, filter=NULL) {

		if (is.null(filter)) {
			filter <- SpatVector$new()
		} else {
			filter <- filter@pntr
		}
		if (is.null(extent)) {
			extent <- double()
		} else {
			extent <- as.vector(ext(extent))
		}
	
		v <- methods::new("SpatVectorCollection")
		v@pntr <- SpatVectorCollection$new(x, layer, query, dialect, extent, filter)	
		v
	}
)



setMethod("svc", signature(x="SpatVector"),
	function(x, ...) {
		r <- methods::new("SpatVectorCollection")
		r@pntr <- SpatVectorCollection$new()
		r@pntr$push_back(x@pntr)
		dots <- list(...)
		if (length(dots) > 0) {
			for (i in 1:length(dots)) {
				if (inherits(dots[[i]], "SpatVector")) {
					r@pntr$push_back(dots[[i]]@pntr)
				} else {
					warn("svc", "cannot add objects of class: ", class(dots[[i]]))
				}
			}
		}
		messages(r, "svc")
	}
)


setMethod("svc", signature(x="sf"),
	function(x) {
		.svc_from_sf(x)
	}
)


setMethod("svc", signature(x="list"),
	function(x) {
		r <- methods::new("SpatVectorCollection")
		r@pntr <- SpatVectorCollection$new()
		for (i in seq_along(x)) {
			if (inherits(x[[i]], "SpatVector")) {
				r@pntr$push_back(x[[i]]@pntr)
			}
		}
		r <- messages(r, "svc")
		names(r) <- names(x)
		r
	}
)


setReplaceMethod("[", c("SpatVectorCollection", "numeric", "missing"),
	function(x, i, j, value) {
		stopifnot(inherits(value, "SpatVector"))
		if (any(!is.finite(i)) || any(i<1)) {
			error("`[<-`", "invalid index")
		}
		i <- sort(i)
		for (j in i) {
			if (j == (length(x)+1)) {
				x@pntr$push_back(value@pntr)
			} else {
				x@pntr$replace(value@pntr, j-1)
			}
		}
		messages(x, "`[<-`")
	}
)


setMethod("[", c("SpatVectorCollection", "numeric", "missing"),
function(x, i, j, drop=TRUE) {
	if (i < 0) {i <- (1:length(x))[i]}
	if (drop && (length(i) == 1)) {
		tptr <- x@pntr$get(i-1)
		x <- methods::new("SpatVector")
		x@pntr <- tptr
	} else {
		x@pntr <- x@pntr$subset(i-1)
	}
	messages(x, "`[`")
})



setMethod("[[", c("SpatVectorCollection", "ANY", "missing"),
function(x, i, j, drop=TRUE) {
	if (inherits(i, "character")) {
		i <- na.omit(match(i, names(x)))
		if (length(i) == 0) {
			error("[[", "no matching names")
		}
	}
	x[i, drop=drop]
})



setMethod("$", c("SpatVectorCollection"),
function(x, name) {
	i <- na.omit(grep(name, names(x)))
	if (length(i) == 0) {
		error("$", "no matching names")
	}
	x[i,drop=TRUE]
})


setMethod("[[", c("SpatVectorCollection", "numeric", "missing"),
function(x, i, j, drop=TRUE) {
	x[i,drop=drop]
})


setMethod("c", signature(x="SpatVector"),
	function(x, ...) {
		svc(x, ...)
	}
)

setMethod("c", signature(x="SpatVectorCollection"),
	function(x, ...) {

		x@pntr <- x@pntr$subset(0:(x@pntr$size()-1) ) ## deep copy
		dots <- list(...)
		for (i in seq_along(dots)) {
			if (inherits(dots[[i]], "SpatVectorCollection")) {
				for (j in 1:length(dots[[i]])) {
					x@pntr$push_back(dots[[i]][[j]]@pntr)
				}
			} else if (inherits(dots[[i]], "SpatVector")) {
				x@pntr$push_back(dots[[i]]@pntr)
			} else {
				error("c", "arguments must be SpatVector or SpatVectorCollection")
			}
		}
		messages(x, "c")
	}
)
