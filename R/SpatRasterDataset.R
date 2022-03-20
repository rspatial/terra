
setMethod("length", signature(x="SpatRasterDataset"),
	function(x) {
		x@ptr$nsds()
	}
)


setMethod("sds", signature(x="character"),
	function(x, ids=0) {

		if (length(x) > 1) {
			r <- lapply(x, rast)
			s <- sds(r)
			names(s) <- tools::file_path_sans_ext(basename(x))
			return(s)
		}

		x <- trimws(x[1])
		if (nchar(x) == 0) {
			error("sds", "provide valid file name(s)")
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
		messages(r, "sds")
	}
)

setMethod("sds", signature(x="missing"),
	function(x) {
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
		messages(r, "sds")
	}
)

setMethod("sds", signature(x="list"),
	function(x) {
		r <- methods::new("SpatRasterDataset")
		r@ptr <- SpatRasterStack$new()
		nms <- names(x)
		if (is.null(nms)) nms <- rep("", length(x))
		for (i in seq_along(x)) {
			if (inherits(x[[i]], "SpatRaster")) {
				r@ptr$add(x[[i]]@ptr, nms[i], "", "", FALSE)
			}
		}
		messages(r, "sds")
	}
)


setMethod("sds", signature(x="array"),
	function(x, crs="", extent=NULL) {
		dims <- dim(x)
		if (length(dims) <= 3) {
			return(sds(rast(x, crs=crs, extent=extent)))
		}
		if (length(dims) > 4) {
			if (length(dims) == 5) {
				if (dims[5] == 1) {
					x <- x[,,,,1]
				} else {
					error("sds,array", "cannot handle an array with 5 dimensions")
				}
			} else {
				error("sds,array", "cannot handle an array with more than 4 dimensions")
			}
		}
		r <- lapply(1:dims[4], function(i) rast(x[,,,i], crs=crs, extent=extent))
		sds(r)
	}
)



setMethod("sds", signature(x="stars"),
	function(x) {
		s <- from_stars(x) 
		if (inherits(s, "SpatRaster")) {
			sds(s)
		} else {
			s
		}
	}
)

setMethod("sds", signature(x="stars_proxy"),
	function(x) {
		s <- from_stars(x) 
		if (inherits(s, "SpatRaster")) {
			sds(s)
		} else {
			s
		}
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
				for (j in 1:(length(dots[[i]]))) {
					if (!x@ptr$add(dots[[i]][[j]]@ptr, sdsnms[j], "", "", FALSE)) {
						messages(x, "c")
					}
				}

			} else if (inherits(dots[[i]], "SpatRaster")) {
				if (is.null(nms)) error("c", "arguments must be named")
				if (!x@ptr$add(dots[[i]]@ptr, nms[i], "", "", FALSE)) {
					messages(x, "c")
				}
			} else {
				error("c", "arguments must be SpatRaster or SpatRasterDataset")
			} 
		}
		messages(x, "c")
	}
)


setReplaceMethod("[", c("SpatRasterDataset","numeric","missing"),
	function(x, i, j, value) {
		if (any(!is.finite(i)) | any(i<1)) {
			error("`[`", "invalid index")
		}
		stopifnot(inherits(value, "SpatRaster"))
		i <- sort(i)
		for (j in i) {
			if (j == (length(x)+1)) {
				x@ptr$add(value@ptr, "", "", "", FALSE)
			} else {
				x@ptr$replace(j-1, value@ptr)
			}
		}
		messages(x, "`[`")
	}
)


setMethod("[", c("SpatRasterDataset", "numeric", "missing"),
function(x, i, j, ... ,drop=TRUE) {
	i <- positive_indices(i, length(x), " [ ")

	if (drop && (length(i) == 1)) {
		ptr <- x@ptr$getsds(i-1)
		x <- rast()
		x@ptr <- ptr
	} else {
		x@ptr <- x@ptr$subset(i-1)
	}
	messages(x, "`[`")
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
	messages(x, "`[`")
})


setMethod("[", c("SpatRasterDataset", "logical", "missing"),
function(x, i, j, ... ,drop=TRUE) {
	x[which(i), ..., drop=drop]
})

setMethod("[", c("SpatRasterDataset", "character", "missing"),
function(x, i, j, ... ,drop=TRUE) {
	i <- match(i, names(x))
	if (any(is.na(i))) {
		error("`[`", "unknown name(s) provided")
	}
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


setMethod("sprc", signature(x="missing"),
	function(x) {
		r <- methods::new("SpatRasterCollection")
		r@ptr <- SpatRasterCollection$new()
		r
	}
)


setMethod("sprc", signature(x="SpatRaster"), 
	function(x, ...) {
		sprc(list(x, ...))
	}
)


setMethod("sprc", signature(x="list"), 
	function(x) {
		n <- length(x)
		ptr <- SpatRasterCollection$new()
		if (n > 0) {
			for (i in 1:n) {
				if (inherits(x[[i]], "SpatRaster")) {
					ptr$add(x[[i]]@ptr)
				} else {
					name <- names(x[[i]])
					cls <- class(x[[i]])
					error("sprc", "list elements should be 'SpatRaster'\n", name, "is of class: ", cls)
				}
			}
		}
		x <- new("SpatRasterCollection")
		x@ptr <- ptr
		x
	}
)


setMethod("length", signature(x="SpatRasterCollection"),
	function(x) {
		x@ptr$length()
	}
)

setMethod("[", c("SpatRasterCollection", "numeric", "missing"),
function(x, i, j, ... ,drop=TRUE) {
	i <- positive_indices(i, length(x), " [ ")
	if (drop && (length(i) == 1)) {
		ptr <- x@ptr$x[i][[1]]
		x <- rast()
		x@ptr <- ptr
	} else {
		s <- x@ptr$x[i]
		ptr <- SpatRasterCollection$new()
		for (i in 1:length(s)) {
			ptr$add(s[[i]])
		}
		x@ptr <- ptr
	}
	messages(x, "`[`")
})
