# Authors: Robert J. Hijmans
# Date : October 2018
# Version 1.0
# License GPL v3


.summarize <- function(x, ..., fun, na.rm=FALSE) {
	dots <- list(...)
	add <- NULL	
	if (length(dots) > 0) {
		cls <- sapply(dots, function(i) inherits(i, "SpatRaster"))
		if (any(cls)) {
			y <- c(dots[cls], x)
			x <- do.call(c, y)
		}
		if (!all(cls)) {
			dots <- dots[!cls]
			i <- sapply(dots, function(x) class(x) %in% c("logical", "integer", "numeric"))
			add <- unlist(dots[i], use.names = FALSE)
		}
	}
		
	if (is.null(add)) {
		x@ptr <- x@ptr$summary(fun, na.rm, .terra_environment$options@ptr)
	} else {
		x@ptr <- x@ptr$summary_numb(fun, add, na.rm, .terra_environment$options@ptr)			
	}
	show_messages(x, fun)
	x		
}

setMethod("Summary", signature(x="SpatRaster"),
	function(x, ..., na.rm=FALSE){
		fun <- as.character(sys.call()[[1L]])
		.summarize(x, ..., fun=fun, na.rm=na.rm)
	}
)


setMethod("mean", signature(x="SpatRaster"),
	function(x, ..., trim=NA, na.rm=FALSE){
		if (!is.na(trim)) {	warning("argument 'trim' is ignored") }
		.summarize(x, ..., fun="mean", na.rm=na.rm)
	}
)


