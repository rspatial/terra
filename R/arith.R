# Author: Robert J. Hijmans
# Date :  September 2018
# Version 1.0
# License GPL v3


setMethod("Arith", signature(e1="SpatRaster", e2="SpatRaster"),
    function(e1, e2){ 
		oper <- as.vector(.Generic)[1]
		stopifnot(oper %in% c("+", "-", "^", "*", "/", "%%")) 
		oper <- ifelse(oper == "%%", "%", oper)
		e1@ptr <- e1@ptr$arith_rast(e2@ptr, oper, .terra_environment$options@ptr)
		show_messages(e1, oper)
	}	
)


setMethod("Arith", signature(e1="SpatRaster", e2="numeric"),
    function(e1, e2){ 
		oper <- as.vector(.Generic)[1]
		stopifnot(oper %in% c("+", "-", "^", "*", "/", "%%")) 
		oper <- ifelse(oper == "%%", "%", oper)
		e1@ptr <- e1@ptr$arith_numb(e2, oper, FALSE, .terra_environment$options@ptr)
		show_messages(e1, oper)				
	}	
)


setMethod("Arith", signature(e1="SpatRaster", e2="missing"),
    function(e1, e2){ 
		methods::callGeneric(0, e1)
	}
)

setMethod("Arith", signature(e1="numeric", e2="SpatRaster"),
    function(e1, e2){ 
		oper <- as.vector(.Generic)[1]
		stopifnot(oper %in% c("+", "-", "^", "*", "/", "%%")) 
		oper <- ifelse(oper == "%%", "%", oper)
		e2@ptr <- e2@ptr$arith_numb(e1, oper, TRUE, .terra_environment$options@ptr)
		show_messages(e2, oper)				
	}	
)


setMethod("Compare", signature(e1="SpatRaster", e2="SpatRaster"),
    function(e1, e2){ 
		oper <- as.vector(.Generic)[1]
		e1@ptr <- e1@ptr$arith_rast(e2@ptr, oper, .terra_environment$options@ptr)
		show_messages(e1, oper)
	}	
)


setMethod("Compare", signature(e1="SpatRaster", e2="numeric"),
    function(e1, e2){ 
		oper <- as.vector(.Generic)[1]
		e1@ptr <- e1@ptr$arith_numb(e2, oper, FALSE, .terra_environment$options@ptr)
		show_messages(e1, oper)
	}
)


setMethod("Compare", signature(e1="numeric", e2="SpatRaster"),
    function(e1, e2){ 
		oper <- as.vector(.Generic)[1]
		e2@ptr <- e2@ptr$arith_numb(e1, oper, TRUE, .terra_environment$options@ptr)
		show_messages(e2, oper)
	}	
)


setMethod("Logic", signature(e1="SpatRaster", e2="SpatRaster"),
    function(e1, e2){ 
		oper <- as.vector(.Generic)[1]
		e1@ptr <- e1@ptr$logic_rast(e2@ptr, oper, .terra_environment$options@ptr)
		show_messages(e1, oper)
	}	
)

setMethod("Logic", signature(e1="SpatRaster", e2="numeric"),
    function(e1, e2){ 
		oper <- as.vector(.Generic)[1]
		e2 <- as.logical(e2)
		e1@ptr <- e1@ptr$logic_numb(e2, oper, .terra_environment$options@ptr)
		show_messages(e1, oper)
	}
)


setMethod("Logic", signature(e1="numeric", e2="SpatRaster"),
    function(e1, e2){ 
		oper <- as.vector(.Generic)[1]
		e1 <- as.logical(e1)
		e2@ptr <- e2@ptr$logic_numb(e1, oper, .terra_environment$options@ptr)
		show_messages(e2, oper)
	}	
)

setMethod("Logic", signature(e1="SpatRaster", e2="logical"),
    function(e1, e2){ 
		oper <- as.vector(.Generic)[1]
		e1@ptr <- e1@ptr$logic_numb(e2, oper, .terra_environment$options@ptr)
		show_messages(e1, oper)
	}
)


setMethod("Logic", signature(e1="logical", e2="SpatRaster"),
    function(e1, e2){ 
		oper <- as.vector(.Generic)[1]
		e2@ptr <- e2@ptr$logic_numb(e1, oper, .terra_environment$options@ptr)
		show_messages(e2, oper)
	}	
)


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


