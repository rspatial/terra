# Author: Robert J. Hijmans
# Date :  September 2018
# Version 1.0
# License GPL v3


setMethod("Arith", signature(e1="SpatExtent", e2="numeric"),
    function(e1, e2){ 
		oper <- as.vector(.Generic)[1]
		if (oper == "%%") { 
			e1@ptr <- e1@ptr$align(e2[1], "")
		} else if (oper == "+") {
			e2 <- rep_len(e2, 4)
			e2[c(1,3)] <- -e2[c(1,3)]
			e1 <- ext(as.vector(e1) + e2)
		} else if (oper == "-") {
			e2 <- rep_len(e2, 4)
			e2[c(1,3)] <- -e2[c(1,3)]
			e1 <- ext(as.vector(e1) - e2)
		} else if (oper == "*") {
			e2 <- abs(rep_len(e2, 4))
			e1 <- as.vector(e1)
			dx <- (e1[2] - e1[1])
			dy <- (e1[4] - e1[3]) 
			mx <- e1[1] + dx/2
			my <- e1[3] + dy/2
			e1[1] <- mx - (dx/2)*e2[1]	
			e1[2] <- mx + (dx/2)*e2[2]		
			e1[3] <- my - (dy/2)*e2[3]	
			e1[4] <- my + (dy/2)*e2[4]	
			e1 <- ext(e1)
		} else if (oper == "/") {
			e2 <- abs(rep_len(e2, 4))
			e1 <- as.vector(e1)
			dx <- (e1[2] - e1[1])
			dy <- (e1[4] - e1[3]) 
			mx <- e1[1] + dx/2
			my <- e1[3] + dy/2
			e1[1] <- mx - dx/(2*e2[1])	
			e1[2] <- mx + dx/(2*e2[2])	
			e1[3] <- my - dy/(2*e2[3])	
			e1[4] <- my + dy/(2*e2[4])	
			e1 <- ext(e1)
		} else {
			stop("only +, - and %% are supported")
		}
		if (!e1@ptr$valid) {
			stop("this would create an invalid extent")
		}
		e1
	}	
)



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


setMethod("!", signature(x="SpatRaster"),
	function(x) {
		x == 0
	}
)	

setMethod("isTRUE", signature(x="SpatRaster"),
	function(x) {
		x <- x == 1
		classify(x, cbind(NA, 0))
	}
)	


setMethod("isFALSE", signature(x="SpatRaster"),
	function(x) {
		x <- x != 1
		classify(x, cbind(NA, 0))
	}
)

setMethod("as.logical", signature(x="SpatRaster"),
	function(x) {
		isTRUE(x)
	}
)	


setMethod("is.na", signature(x="SpatRaster"),
	function(x) {
		x@ptr <- x@ptr$isnan(.terra_environment$options@ptr)
		show_messages(x, "is.na")
	}
)	


setMethod("is.nan", signature(x="SpatRaster"),
	function(x) {
		x@ptr <- x@ptr$isnan(.terra_environment$options@ptr)
		show_messages(x, "is.nan")
	}
)	


setMethod("is.finite", signature(x="SpatRaster"),
	function(x) {
		x@ptr <- x@ptr$isfinite(.terra_environment$options@ptr)
		show_messages(x, "is.finite")
	}
)	

setMethod("is.infinite", signature(x="SpatRaster"),
	function(x) {
		x@ptr <- x@ptr$isinfinite(.terra_environment$options@ptr)
		show_messages(x, "is.infinite")
	}
)	


.summarize <- function(x, ..., fun, na.rm=FALSE) {
	dots <- list(...)
	add <- NULL	
	cls <- FALSE
	if (length(dots) > 0) {
		cls <- sapply(dots, function(i) inherits(i, "SpatRaster"))
		if (!all(cls)) {
			dots <- dots[!cls]
			i <- sapply(dots, function(x) class(x) %in% c("logical", "integer", "numeric"))
			add <- unlist(dots[i], use.names = FALSE)
		}
	}
	if (any(cls)) {
		x <- sds(c(list(x), dots[cls]))
	} 

	r <- rast()
	if (is.null(add)) {
		r@ptr <- x@ptr$summary(fun, na.rm, .terra_environment$options@ptr)
	} else {
		r@ptr <- x@ptr$summary_numb(fun, add, na.rm, .terra_environment$options@ptr)			
	}
	show_messages(r, fun)
	r
}


setMethod("which.max", "SpatRaster",  
	function(x) { 
		x@ptr <- x@ptr$summary("which.max", TRUE, .terra_environment$options@ptr)
		show_messages(x, "which.max")
	}
)

setMethod("which.min", "SpatRaster",  
	function(x) { 
		x@ptr <- x@ptr$summary("which.min", TRUE, .terra_environment$options@ptr)
		show_messages(x, "which.min")
	}
)



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

setMethod("median", signature(x="SpatRaster"),
	function(x, na.rm=FALSE, ...){
		.summarize(x, ..., fun="median", na.rm=na.rm)
	}
)


setMethod("Compare", signature(e1="SpatExtent", e2="SpatExtent"),
    function(e1, e2){ 
		oper <- as.vector(.Generic)[1]
		if (!(oper %in% c("==", "!=", ">", "<", ">=", "<="))) {
			stop(paste(oper, "not implemented for SpatExtent"))
		}
		return( e1@ptr$compare(e2@ptr, oper, 0.000001) )
	}	
)


setMethod("stdev", signature(x="SpatRaster"),
	function(x, ..., na.rm=FALSE){
		.summarize(x, ..., fun="stdev", na.rm=na.rm)
	}
)



setMethod("modal", signature("SpatRaster"), 
	function(x, ..., ties="first", na.rm=FALSE, filename="", overwrite=FALSE, wopt=list()) { 
		opt <- .runOptions(filename, overwrite,wopt)
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
			add <- c(.5)[0]
		}
		x@ptr <- x@ptr$modal(add, ties[1], na.rm[1], opt)
		show_messages(x, "modal")		
	}
)

