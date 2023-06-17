
setMethod("length", signature(x="SpatVectorCollection"),
	function(x) {
		x@ptr$size()
	}
)

setMethod("svc", signature(x="missing"),
	function(x) {
		v <- methods::new("SpatVectorCollection")
		v@ptr <- SpatVectorCollection$new()
		v
	}
)

setMethod("svc", signature(x="character"),
	function(x, layer="", query="", extent=NULL, filter=NULL) {

		if (is.null(filter)) {
			filter <- SpatVector$new()
		} else {
			filter <- filter@ptr
		}
		if (is.null(extent)) {
			extent <- double()
		} else {
			extent <- as.vector(ext(extent))
		}
	
		v <- methods::new("SpatVectorCollection")
		v@ptr <- SpatVectorCollection$new(x, layer, query, extent, filter)	
		v
	}
)



setMethod("svc", signature(x="SpatVector"),
	function(x, ...) {
		r <- methods::new("SpatVectorCollection")
		r@ptr <- SpatVectorCollection$new()
		r@ptr$push_back(x@ptr)
		dots <- list(...)
		if (length(dots) > 0) {
			for (i in 1:length(dots)) {
				if (inherits(dots[[i]], "SpatVector")) {
					r@ptr$push_back(dots[[i]]@ptr)
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
		r@ptr <- SpatVectorCollection$new()
		for (i in seq_along(x)) {
			if (inherits(x[[i]], "SpatVector")) {
				r@ptr$push_back(x[[i]]@ptr)
			}
		}
		messages(r, "svc")
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
				x@ptr$push_back(value@ptr)
			} else {
				x@ptr$replace(value@ptr, j-1)
			}
		}
		messages(x, "`[<-`")
	}
)


setMethod("[", c("SpatVectorCollection", "numeric", "missing"),
function(x, i, j, drop=TRUE) {
	if (i < 0) {i <- (1:length(x))[i]}
	if (drop && (length(i) == 1)) {
		ptr <- x@ptr$get(i-1)
		x <- methods::new("SpatVector")
		x@ptr <- ptr
	} else {
		x@ptr <- x@ptr$subset(i-1)
	}
	messages(x, "`[`")
})

setMethod("[[", c("SpatVectorCollection", "numeric", "missing"),
function(x, i, drop=TRUE) {
	x[i,drop=drop]
})


setMethod("c", signature(x="SpatVector"),
	function(x, ...) {
		svc(x, ...)
	}
)

setMethod("c", signature(x="SpatVectorCollection"),
	function(x, ...) {

		x@ptr <- x@ptr$subset(0:(x@ptr$size()-1) ) ## deep copy
		dots <- list(...)
		for (i in seq_along(dots)) {
			if (inherits(dots[[i]], "SpatVectorCollection")) {
				for (j in 1:length(dots[[i]])) {
					x@ptr$push_back(dots[[i]][[j]]@ptr)
				}
			} else if (inherits(dots[[i]], "SpatVector")) {
				x@ptr$push_back(dots[[i]]@ptr)
			} else {
				error("c", "arguments must be SpatVector or SpatVectorCollection")
			}
		}
		messages(x, "c")
	}
)
