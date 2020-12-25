
setMethod("relate", signature(x="SpatVector", y="SpatVector"), 
	function(x, y, relation) {
		relation <- tolower(relation)
		relation <- match.arg(relation, c("intersects", "touches", "crosses", "overlaps", "within", "contains", "covers", "coveredby", "disjoint"))
		out <- x@ptr$relate_between(y@ptr, relation)
		x <- messages(x, "relate")
		out[out == 2] <- NA
		matrix(as.logical(out), nrow=nrow(x), byrow=TRUE)
	}
)


setMethod("relate", signature(x="SpatVector", y="missing"), 
	function(x, y, relation) {
		relation <- tolower(relation)
		relation <- match.arg(relation, c("intersects", "touches", "crosses", "overlaps", "within", "contains", "covers", "coveredby", "disjoint"))
		out <- x@ptr$relate_within(relation)
		x <- messages(x, "relate")
		out[out == 2] <- NA
		matrix(as.logical(out), nrow=nrow(x), byrow=TRUE)
	}
)

setMethod("relate", signature(x="SpatVector", y="ANY"), 
	function(x, y, relation) {
		yy <- try(vect(y), silent=TRUE)
		if (!inherits(yy, "SpatVector")) {
			yy <- try(ext(y), silent=TRUE)
			if (!inherits(yy, "SpatExtent")) {
				stop("cannot use argument 'y'")
			}
		}
		relate(x, yy, relation)
	}
)

setMethod("relate", signature(x="ANY", y="SpatVector"), 
	function(x, y, relation) {
		xx <- try(vect(x), silent=TRUE)
		if (!inherits(xx, "SpatVector")) {
			xx <- try(ext(x), silent=TRUE)
			if (!inherits(xx, "SpatExtent")) {
				stop("cannot use argument 'x'")
			}
		}
		relate(xx, y, relation)
	}
)

setMethod("relate", signature(x="ANY", y="ANY"), 
	function(x, y, relation) {
		xx <- try(vect(x), silent=TRUE)
		if (!inherits(xx, "SpatVector")) {
			xx <- try(ext(x), silent=TRUE)
			if (!inherits(xx, "SpatExtent")) {
				stop("cannot use argument 'x'")
			}
		}
		yy <- try(vect(y), silent=TRUE)
		if (!inherits(yy, "SpatVector")) {
			yy <- try(ext(y), silent=TRUE)
			if (!inherits(yy, "SpatExtent")) {
				stop("cannot use argument 'y'")
			}
		}
		relate(xx, yy, relation)
	}
)
