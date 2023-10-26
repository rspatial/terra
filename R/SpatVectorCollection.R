
setMethod("length", signature(x="SpatVectorCollection"),
	function(x) {
		x@cpp$size()
	}
)

setMethod("svc", signature(x="missing"),
	function(x) {
		v <- methods::new("SpatVectorCollection")
		v@cpp <- SpatVectorCollection$new()
		v
	}
)

setMethod("svc", signature(x="character"),
	function(x, layer="", query="", extent=NULL, filter=NULL) {

		if (is.null(filter)) {
			filter <- SpatVector$new()
		} else {
			filter <- filter@cpp
		}
		if (is.null(extent)) {
			extent <- double()
		} else {
			extent <- as.vector(ext(extent))
		}
	
		v <- methods::new("SpatVectorCollection")
		v@cpp <- SpatVectorCollection$new(x, layer, query, extent, filter)	
		v
	}
)



setMethod("svc", signature(x="SpatVector"),
	function(x, ...) {
		r <- methods::new("SpatVectorCollection")
		r@cpp <- SpatVectorCollection$new()
		r@cpp$push_back(x@cpp)
		dots <- list(...)
		if (length(dots) > 0) {
			for (i in 1:length(dots)) {
				if (inherits(dots[[i]], "SpatVector")) {
					r@cpp$push_back(dots[[i]]@cpp)
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
		r@cpp <- SpatVectorCollection$new()
		for (i in seq_along(x)) {
			if (inherits(x[[i]], "SpatVector")) {
				r@cpp$push_back(x[[i]]@cpp)
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
				x@cpp$push_back(value@cpp)
			} else {
				x@cpp$replace(value@cpp, j-1)
			}
		}
		messages(x, "`[<-`")
	}
)


setMethod("[", c("SpatVectorCollection", "numeric", "missing"),
function(x, i, j, drop=TRUE) {
	if (i < 0) {i <- (1:length(x))[i]}
	if (drop && (length(i) == 1)) {
		ptr <- x@cpp$get(i-1)
		x <- methods::new("SpatVector")
		x@cpp <- ptr
	} else {
		x@cpp <- x@cpp$subset(i-1)
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

		x@cpp <- x@cpp$subset(0:(x@cpp$size()-1) ) ## deep copy
		dots <- list(...)
		for (i in seq_along(dots)) {
			if (inherits(dots[[i]], "SpatVectorCollection")) {
				for (j in 1:length(dots[[i]])) {
					x@cpp$push_back(dots[[i]][[j]]@cpp)
				}
			} else if (inherits(dots[[i]], "SpatVector")) {
				x@cpp$push_back(dots[[i]]@cpp)
			} else {
				error("c", "arguments must be SpatVector or SpatVectorCollection")
			}
		}
		messages(x, "c")
	}
)
