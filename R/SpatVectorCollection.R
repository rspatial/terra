
setMethod("length", signature(x="SpatVectorCollection"),
	function(x) {
		x@pnt$size()
	}
)

setMethod("svc", signature(x="missing"),
	function(x) {
		v <- methods::new("SpatVectorCollection")
		v@pnt <- SpatVectorCollection$new()
		v
	}
)

setMethod("svc", signature(x="character"),
	function(x, layer="", query="", extent=NULL, filter=NULL, crs="") {

		if (is.null(filter)) {
			filter <- SpatVector$new()
		} else {
			if (proxy) {
				error("vect", "you cannot use 'filter' when proxy=TRUE")
			}
			filter <- filter@pnt
		}
		if (is.null(extent)) {
			extent <- double()
		} else {
			extent <- as.vector(ext(extent))
		}
	
		v <- methods::new("SpatVectorCollection")
		v@pnt <- SpatVectorCollection$new(x, layer, query, extent, filter)
		
		if (isTRUE(crs != "")) {
			crs(v, warn=FALSE) <- crs
		}

		v
	}
)



setMethod("svc", signature(x="SpatVector"),
	function(x, ...) {
		r <- methods::new("SpatVectorCollection")
		r@pnt <- SpatVectorCollection$new()
		r@pnt$push_back(x@pnt)
		dots <- list(...)
		if (length(dots) > 0) {
			for (i in 1:length(dots)) {
				if (inherits(dots[[i]], "SpatVector")) {
					r@pnt$push_back(dots[[i]]@pnt)
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
		r@pnt <- SpatVectorCollection$new()
		for (i in seq_along(x)) {
			if (inherits(x[[i]], "SpatVector")) {
				r@pnt$push_back(x[[i]]@pnt)
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
				x@pnt$push_back(value@pnt)
			} else {
				x@pnt$replace(value@pnt, j-1)
			}
		}
		messages(x, "`[<-`")
	}
)


setMethod("[", c("SpatVectorCollection", "numeric", "missing"),
function(x, i, j, drop=TRUE) {
	if (i < 0) {i <- (1:length(x))[i]}
	if (drop && (length(i) == 1)) {
		ptr <- x@pnt$get(i-1)
		x <- methods::new("SpatVector")
		x@pnt <- ptr
	} else {
		x@pnt <- x@pnt$subset(i-1)
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

		x@pnt <- x@pnt$subset(0:(x@pnt$size()-1) ) ## deep copy
		dots <- list(...)
		for (i in seq_along(dots)) {
			if (inherits(dots[[i]], "SpatVectorCollection")) {
				for (j in 1:length(dots[[i]])) {
					x@pnt$push_back(dots[[i]][[j]]@pnt)
				}
			} else if (inherits(dots[[i]], "SpatVector")) {
				x@pnt$push_back(dots[[i]]@pnt)
			} else {
				error("c", "arguments must be SpatVector or SpatVectorCollection")
			}
		}
		messages(x, "c")
	}
)
