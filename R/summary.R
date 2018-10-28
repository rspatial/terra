# Authors: Robert J. Hijmans
# Date : October 2018
# Version 1.0
# Licence GPL v3

setMethod("Summary", signature(x='SpatRaster'),
	function(x, ..., na.rm=FALSE){
		dots <- list(...)
		add <- NULL	
		if (length(dots) > 0) {
			cls <- sapply(y, function(i) inherits(i, "SpatRaster"))
			if (any(cls)) {
				y <- c(dots[cls], x)
				x <- do.call(c, y)
			}
			if (!all(cls)) {
				dots <- dots[!cls]
				i <- sapply(dots, function(x) class(x) %in% c('logical', 'integer', 'numeric'))
				add <- unlist(dots[i], use.names = FALSE)
			}
		}
		
		fun <- as.character(sys.call()[[1L]])
		if (fun == "sum") {
			if (is.null(add)) {
				x@ptr <- x@ptr$summary("sum", na.rm, "", TRUE)
			} else {
				x@ptr <- x@ptr$summary_numb("sum", add, na.rm, "", TRUE)			
			}
		} else if (fun == 'min') {
			if (is.null(add)) {
				x@ptr <- x@ptr$summary("min", na.rm, "", TRUE)
			} else {
				x@ptr <- x@ptr$summary_numb("min", add, na.rm, "", TRUE)			
			}
		} else if (fun == 'max') {
			if (is.null(add)) {
				x@ptr <- x@ptr$summary("max", na.rm, "", TRUE)
			} else {
				x@ptr <- x@ptr$summary_numb("max", add, na.rm, "", TRUE)			
			}
		} else if (fun == 'range') {
			if (is.null(add)) {
				x@ptr <- x@ptr$summary("range", na.rm, "", TRUE)
			} else {
				x@ptr <- x@ptr$summary_numb("range", add, na.rm, "", TRUE)			
			}
		}
		.messages(x, fun)
		x		
	}
)


setMethod("mean", signature(x='SpatRaster'),
	function(x, ..., trim=NA, na.rm=FALSE){

		if (!is.na(trim)) {	warning("argument 'trim' is ignored") }

		add <- NULL	
		dots <- list(...)
		if (length(dots) > 0) {
			cls <- sapply(y, function(i) inherits(i, "SpatRaster"))
			if (any(cls)) {
				y <- c(dots[cls], x)
				x <- do.call(c, y)
			}
			if (!all(cls)) {
				dots <- dots[!cls]
				i <- sapply(dots, function(x) class(x) %in% c('logical', 'integer', 'numeric'))
				add <- unlist(dots[i], use.names = FALSE)
			}
		}
		
		if (is.null(add)) {
			x@ptr <- x@ptr$summary("mean", na.rm, "", TRUE)
		} else {
			x@ptr <- x@ptr$summary_numb("mean", add, na.rm, "", TRUE)			
		}
		.messages(x, "mean")
		x		
	}
)


