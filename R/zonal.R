
setMethod("zonal", signature(x="SpatRaster", z="SpatRaster"), 
	function(x, z, fun="mean", ...)  {
		txtfun <- .makeTextFun(match.fun(fun))
		if (inherits(txtfun, "character")) { 
			if (txtfun %in% c("max", "min", "mean", "sum")) {
				na.rm <- isTRUE(list(...)$na.rm)
				opt <- .getOptions()
				ptr <- x@ptr$zonal(z@ptr, txtfun, na.rm, opt)
				show_messages(ptr, "zonal")
				return( .getSpatDF(ptr) )
			}		
		} 
		
		#else 
		nl <- nlyr(x)
		res <- list()
		z <- values(z)
		nms <- names(x)
		for (i in 1:nl) {
			d <- stats::aggregate(values(x[[i]]), list(zone=z), fun, ...)
			colnames(d)[2] <- nms[i]
			res[[i]] <- d
		}
		r <- res[[1]]
		if (nl > 1) {
			for (i in 2:nl) {
				r <- merge(r, res[[i]])
			}
		}
		r
	}
)


setMethod("global", signature(x="SpatRaster"), 
	function(x, fun="mean", weights=NULL, ...)  {

		nms <- names(x)
		nms <- make.unique(nms)
		txtfun <- .makeTextFun(fun)

		#if (txtfun == '.Primitive(\"range\")') txtfun <- "range" 

		opt <- .getOptions()
		if (!is.null(weights)) {
			stopifnot(inherits(weights, "SpatRaster"))
			stopifnot(txtfun %in% c("mean", "sum"))
			na.rm <- isTRUE(list(...)$na.rm)
			ptr <- x@ptr$global_weighted_mean(weights@ptr, txtfun, na.rm, opt)
			show_messages(ptr, "global")
			res <- (.getSpatDF(ptr))
			rownames(res) <- nms
			return(res)
		}
		
		if (inherits(txtfun, "character")) { 
			if (txtfun %in% c("max", "min", "mean", "sum", "range", "rms")) {
				na.rm <- isTRUE(list(...)$na.rm)
				ptr <- x@ptr$global(txtfun, na.rm, opt)
				show_messages(ptr, "global")
				res <- .getSpatDF(ptr)

				rownames(res) <- nms
				return(res)
			}
		}

		nl <- nlyr(x)
		res <- list()
		for (i in 1:nl) {
			res[[i]] <- fun(values(x[[i]]), ...)
		}
		res <- do.call(rbind,res)
		res <- data.frame(res)
		if (ncol(res) > 1) {
			colnames(res) <- paste0("global_", 1:ncol(res))			
		} else {
			colnames(res) <- "global"
		}
		rownames(res) <- nms
		res
	}
)
