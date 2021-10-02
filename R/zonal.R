
setMethod("zonal", signature(x="SpatRaster", z="SpatRaster"), 
	function(x, z, fun="mean", ..., as.raster=FALSE, filename="", wopt=list())  {
		txtfun <- .makeTextFun(match.fun(fun))
		if (nlyr(z) > 1) {
			z <- z[[1]]
		}
		if (inherits(txtfun, "character")) { 
			if (txtfun %in% c("max", "min", "mean", "sum")) {
				na.rm <- isTRUE(list(...)$na.rm)
				opt <- spatOptions()
				ptr <- x@ptr$zonal(z@ptr, txtfun, na.rm, opt)
				messages(ptr, "zonal")
				out <- .getSpatDF(ptr)
			}
		} else {
			nl <- nlyr(x)
			res <- list()
			z <- values(z)
			nms <- names(x)
			for (i in 1:nl) {
				d <- stats::aggregate(values(x[[i]]), list(zone=z), fun, ...)
				colnames(d)[2] <- nms[i]
				res[[i]] <- d
			}
			out <- res[[1]]
			if (nl > 1) {
				for (i in 2:nl) {
					out <- merge(out, res[[i]])
				}
			}
		}

		if (as.raster) {
			subst(z, out[,1], out[,-1], filename=filename, wopt=wopt)
		} else {
			if (is.factor(z)) {
				levs <- levels(z)[[1]]
				out$zone <- levs[out$zone+1]
			}
			colnames(out)[1] <- names(z)
			out
		}
	}
)


setMethod("global", signature(x="SpatRaster"), 
	function(x, fun="mean", weights=NULL, ...)  {

		nms <- names(x)
		nms <- make.unique(nms)
		txtfun <- .makeTextFun(fun)

		opt <- spatOptions()
		if (!is.null(weights)) {
			stopifnot(inherits(weights, "SpatRaster"))
			stopifnot(txtfun %in% c("mean", "sum"))
			na.rm <- isTRUE(list(...)$na.rm)
			ptr <- x@ptr$global_weighted_mean(weights@ptr, txtfun, na.rm, opt)
			messages(ptr, "global")
			res <- (.getSpatDF(ptr))
			rownames(res) <- nms
			return(res)
		}

		if (inherits(txtfun, "character")) { 
			if (txtfun %in% c("max", "min", "mean", "sum", "range", "rms", "sd", "sdpop")) {
				na.rm <- isTRUE(list(...)$na.rm)
				ptr <- x@ptr$global(txtfun, na.rm, opt)
				messages(ptr, "global")
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
