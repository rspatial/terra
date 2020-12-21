
setMethod("length", signature(x="SpatRasterDataset"),
	function(x) {
		x@ptr$size()
	}
)

setMethod("svc", signature(x="missing"),
	function(x, ...) {
		v <- methods::new("SpatVectorCollection")
		v@ptr <- SpatRasterStack$new()
		v
	}
)


setMethod("svc", signature(x="SpatVector"),
	function(x, ...) {
		r <- methods::new("SpatVectorCollection")
		r@ptr <- SpatVectorCollection$new()
		r@ptr$push_back(x@ptr)
		messages(r, "svc")
	}
)

setMethod("sds", signature(x="list"),
	function(x, ...) {
		r <- methods::new("SpatVectorCollection")
		r@ptr <- SpatVectorCollection$new()
		for (i in seq_along(x)) {
			if (inherits(x[[i]], "SpatVector")) {
				r@ptr$push_back(x[[i]]@ptr)
			}
		}	
		messages(r, "sds")
	}
)


setReplaceMethod("[", c("SpatVectorCollection","numeric","missing"),
	function(x, i, j, value) {
		stopifnot(inherits(value, "SpatVector"))
		if (any(!is.finite(i)) | any(i<1)) {
			error(" [,SpatVectorCollection", "invalid index")
		}
		if (length(i) > 1) {
			error(" [,SpatVectorCollection", "you can only replace one sub-dataset at a time")		
		}
		x@ptr$replace(i-1, value@ptr)
		messages(x, "[")
	}
)


setMethod("[", c("SpatVectorCollection", "numeric", "missing"),
function(x, i, j, ... ,drop=TRUE) {
	if (i<0) {i <- (1:length(x))[i]}
	if (drop && (length(i) == 1)) {
		ptr <- x@ptr$get(i-1)
		x <- methods::new("SpatVector")
		x@ptr <- ptr
	} else {
		x@ptr <- x@ptr$subset(i-1)
	}
	messages(x, "[")
})

