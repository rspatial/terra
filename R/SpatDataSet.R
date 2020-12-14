
setMethod("length", signature(x="SpatRasterDataset"),
	function(x) {
		x@ptr$nsds()
	}
)

setMethod("sds", signature(x="character"),
	function(x, ids=0, ...) {
		
		if (length(x) > 1) {
			r <- lapply(x, rast)
			s <- sds(r)	
			names(s) <- tools::file_path_sans_ext(basename(x))
			return(s)
		}
		
		x <- trimws(x[1])
		if (nchar(x) == 0) {
			stop("provide valid file name(s)")
		}
		f <- .fullFilename(x)
		r <- methods::new("SpatRasterDataset")
		ids <- round(ids)-1
		if (ids[1] < 0) {
			useids <- FALSE
		} else {
			useids <- TRUE
		}
		r@ptr <- SpatRasterStack$new(f, ids, useids)
		show_messages(r, "sds")
	}
)

setMethod("sds", signature(x="missing"),
	function(x, ...) {
		r <- methods::new("SpatRasterDataset")
		r@ptr <- SpatRasterStack$new()
		r
	}
)


setMethod("sds", signature(x="SpatRaster"),
	function(x, ...) {
		r <- methods::new("SpatRasterDataset")
		r@ptr <- SpatRasterStack$new()
		r@ptr$add(x@ptr, varnames(x)[1], longnames(x)[1], units(x)[1], FALSE)
		dots <- list(...)
		nms <- names(dots)
		if (is.null(nms)) nms = ""
		nms <- rep_len(nms, length(dots))
		for (i in seq_along(dots)) {
			if (inherits(dots[[i]], "SpatRaster")) {
				r@ptr$add(dots[[i]]@ptr, nms[i], "", "", FALSE)
			}
		}	
		show_messages(r, "sds")
	}
)

setMethod("sds", signature(x="list"),
	function(x, ...) {
		r <- methods::new("SpatRasterDataset")
		r@ptr <- SpatRasterStack$new()
		nms <- names(x)
		if (is.null(nms)) nms <- rep("", length(x))
		for (i in seq_along(x)) {
			if (inherits(x[[i]], "SpatRaster")) {
				r@ptr$add(x[[i]]@ptr, nms[i], "", "", FALSE)
			}
		}	
		show_messages(r, "sds")
	}
)

setMethod("c", signature(x="SpatRasterDataset"), 
	function(x, ...) {
		
		x@ptr <- x@ptr$subset((1:x@ptr$nsds()) -1 ) # why? make a copy?
	 	
		dots <- list(...)
		nms <- names(dots)
		
		for (i in seq_along(dots)) {
			if (inherits(dots[[i]], "SpatRasterDataset")) {
				sdsnms <- names(dots[[i]])
				for (j in 1:x@ptr$nsds()) {
					if (!x@ptr$add(dots[[i]][[j]]@ptr, sdsnms[j], "", "", FALSE)) {
						show_messages(x, "c")		
					}
				}
			
			} else if (inherits(dots[[i]], "SpatRaster")) {
				if (is.null(nms)) stop("arguments must be named")
				if (!x@ptr$add(dots[[i]]@ptr, nms[i], "", "", FALSE)) {
					show_messages(x, "c")		
				}
			} else {
				stop("arguments must be SpatRaster or SpatRasterDataset")
			} 
		}
		show_messages(x, "c")		
	}
)


setReplaceMethod("[", c("SpatRasterDataset","numeric","missing"),
	function(x, i, j, value) {
		if (any(!is.finite(i)) | any(i<1)) {
			stop("invalid index")
		}
		if (length(i) > 1) {
			stop("you can only replace one sub-dataset at a time")		
		}
		stopifnot(inherits(value, "SpatRaster"))
		x@ptr$replace(i-1, value@ptr)
		show_messages(x, "[")
	}
)


setMethod("[", c("SpatRasterDataset", "numeric", "missing"),
function(x, i, j, ... ,drop=TRUE) {
	if (i<0) {i <- (1:length(x))[i]}
	if (drop && (length(i) == 1)) {
		ptr <- x@ptr$getsds(i-1)
		x <- rast()
		x@ptr <- ptr
	} else {
		x@ptr <- x@ptr$subset(i-1)
	}
	show_messages(x, "[")
})

setMethod("[", c("SpatRasterDataset", "numeric", "numeric"),
function(x, i, j, ... ,drop=TRUE) {
	y <- x[i,drop=drop]
	if (inherits(y, "SpatRaster")) {
		return(y[[j]])
	}
	nd <- y@ptr$nsds()
	x@ptr <- SpatRasterStack$new()
	nms <- y@ptr$names
	for (k in seq_along(1:nd)) {
		r <- y[k][[j]]
		x@ptr$add(r@ptr, nms[k], "", "", FALSE)
	}
	show_messages(x, "[")
})


setMethod("[", c("SpatRasterDataset", "logical", "missing"),
function(x, i, j, ... ,drop=TRUE) {
	x[which(i), ..., drop=drop]
})

setMethod("[", c("SpatRasterDataset", "character", "missing"),
function(x, i, j, ... ,drop=TRUE) {
	i <- match(i, names(x))
	if (any(is.na(i))) {stop("unknown name(s) provided")}
	x[i, ..., drop=drop]
})

setMethod("[[", c("SpatRasterDataset", "ANY", "ANY"),
function(x, i, j, ... ,drop=TRUE) {
	mi <- missing(i)
	mj <- missing(j)
	
	if ((mi) && (mj)) {
		`[`(x, ..., drop=drop)
	} else if (mi) {
		`[`(x, j=j, ..., drop=drop)
	} else if (mj) {
		`[`(x, i=i, ..., drop=drop)
	} else {
		`[`(x, i=i, j=j, ..., drop=drop)
	}
})


setMethod("$", "SpatRasterDataset",  
	function(x, name) { 
		x[name] 
	}
)

