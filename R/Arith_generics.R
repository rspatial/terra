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
			error(oper, "only +, -, *, / and %% are supported")
		}
		if (!e1@ptr$valid) {
			error(oper, "this would create an invalid extent")
		}
		e1
	}
)


setMethod("Arith", signature(e1="numeric", e2="SpatExtent"),
    function(e1, e2) {
		oper <- as.vector(.Generic)[1]
		if (oper == "%%") { 
			error("%%", "only 'Spatextent %% numeric' (in that order) is supported")
		} else if (oper == "+") {
			return(e2 + e1)
		} else if (oper == "-") {
			error("-", "only 'Spatextent - numeric' (in that order) is supported")
		} else if (oper == "*") {
			return(e2 * e1)
		} else if (oper == "/") {
			error("/", "only 'Spatextent / numeric' (in that order) is supported")
		} else {
			error(oper, "only +, -, *, / and %% are supported")
		}
		if (!e1@ptr$valid) {
			error(oper, "this would create an invalid extent")
		}
		e1
	}
)


setMethod("Arith", signature(e1="SpatExtent", e2="SpatExtent"),
    function(e1, e2){ 
		oper <- as.vector(.Generic)[1]
		e1 = ext(as.vector(e1)) # deep copy
		if (oper == "+") { 
			e1@ptr$union(e2@ptr)
		} else if (oper == "*") {
			e1@ptr$intersect(e2@ptr)
		} else if (oper == "/") {
			d <- c(diff(e1[1:2]) / diff(e2[1:2]), diff(e1[3:4]) / diff(e2[3:4])) 
			names(d) <- c("x", "y")
			return(d)
		} else {
			error(oper, "only +, *, and / are supported for SpatExtent")
		}
		if (!e1@ptr$valid) {
			error(oper, "this would create an invalid extent")
		}
		e1
	}
)



setMethod("Arith", signature(e1="SpatVector", e2="SpatVector"),
    function(e1, e2){ 
		oper <- as.vector(.Generic)[1]
		if (oper == "+") { 
			e1@ptr <- e1@ptr$union(e2@ptr)
		} else if (oper == "*") {
			e1@ptr <- e1@ptr$intersect(e2@ptr)
		} else if (oper == "-") {
			e1@ptr <- e1@ptr$erase(e2@ptr)
		} else {
			error(oper, "only operators +, *, and - are supported for SpatVector")
		}
		messages(e1, oper)
	}
)




setMethod("Arith", signature(e1="SpatRaster", e2="SpatRaster"),
    function(e1, e2){ 
		opt <- spatOptions()
		oper <- as.vector(.Generic)[1]
		stopifnot(oper %in% c("+", "-", "^", "*", "/", "%%")) 
		oper <- ifelse(oper == "%%", "%", oper)
		e1@ptr <- e1@ptr$arith_rast(e2@ptr, oper, opt)
		messages(e1, oper)
	}
)


setMethod("Arith", signature(e1="SpatRaster", e2="numeric"),
    function(e1, e2){
		oper <- as.vector(.Generic)[1]
		stopifnot(oper %in% c("+", "-", "^", "*", "/", "%%")) 
		opt <- spatOptions()
		oper <- ifelse(oper == "%%", "%", oper)
		e1@ptr <- e1@ptr$arith_numb(e2, oper, FALSE, opt)
		messages(e1, oper)
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
		opt <- spatOptions()
		oper <- ifelse(oper == "%%", "%", oper)
		e2@ptr <- e2@ptr$arith_numb(e1, oper, TRUE, opt)
		messages(e2, oper)
	}
)


setMethod("Compare", signature(e1="SpatRaster", e2="SpatRaster"),
    function(e1, e2){ 
		opt <- spatOptions()
		oper <- as.vector(.Generic)[1]
		e1@ptr <- e1@ptr$arith_rast(e2@ptr, oper, opt)
		messages(e1, oper)
	}
)


setMethod("Compare", signature(e1="SpatRaster", e2="numeric"),
    function(e1, e2){ 
		oper <- as.vector(.Generic)[1]
		opt <- spatOptions()
		e1@ptr <- e1@ptr$arith_numb(e2, oper, FALSE, opt)
		messages(e1, oper)
	}
)


setMethod("Compare", signature(e1="numeric", e2="SpatRaster"),
    function(e1, e2){ 
		oper <- as.vector(.Generic)[1]
		opt <- spatOptions()
		e2@ptr <- e2@ptr$arith_numb(e1, oper, TRUE, opt)
		messages(e2, oper)
	}
)


setMethod("Compare", signature(e1="SpatRaster", e2="character"),
    function(e1, e2){ 
		oper <- as.vector(.Generic)[1]
		if (!is.factor(e1)) {
			error(oper, "SpatRaster is not categorical")		
		}
		if (oper != "==") {
			error(oper, "only '==' is supported with categorical comparisons")
		}
		if (nlyr(e1) != 1) {
			error(oper, "categorical comparisons only supported for single layer SpatRaster")
		}
		if (length(e2) != 1) {
			error(oper, "comparisons only supported for single values (see %in% and match)")
		}
		
		e2 <- match(e2, levels(e1)[[1]])
		if (is.na(e2)) return (e1 * 0)
		opt <- spatOptions()
		e1@ptr <- e1@ptr$arith_numb(e2, oper, TRUE, opt)
		messages(e1, oper)
	}
)



setMethod("Logic", signature(e1="SpatRaster", e2="SpatRaster"),
    function(e1, e2){ 
		oper <- as.vector(.Generic)[1]
		opt <- spatOptions()
		e1@ptr <- e1@ptr$logic_rast(e2@ptr, oper, opt)
		messages(e1, oper)
	}
)

setMethod("Logic", signature(e1="SpatRaster", e2="numeric"),
    function(e1, e2){ 
		opt <- spatOptions()
		oper <- as.vector(.Generic)[1]
		e2 <- as.logical(e2)
		e1@ptr <- e1@ptr$logic_numb(e2, oper, opt)
		messages(e1, oper)
	}
)


setMethod("Logic", signature(e1="numeric", e2="SpatRaster"),
    function(e1, e2){ 
		opt <- spatOptions()
		oper <- as.vector(.Generic)[1]
		e1 <- as.logical(e1)
		e2@ptr <- e2@ptr$logic_numb(e1, oper, opt)
		messages(e2, oper)
	}
)

setMethod("Logic", signature(e1="SpatRaster", e2="logical"),
    function(e1, e2){ 
		opt <- spatOptions()
		oper <- as.vector(.Generic)[1]
		e1@ptr <- e1@ptr$logic_numb(e2, oper, opt)
		messages(e1, oper)
	}
)


setMethod("Logic", signature(e1="logical", e2="SpatRaster"),
    function(e1, e2){ 
		opt <- spatOptions()
		oper <- as.vector(.Generic)[1]
		e2@ptr <- e2@ptr$logic_numb(e1, oper, opt)
		messages(e2, oper)
	}
)


setMethod("!", signature(x="SpatRaster"),
	function(x) {
		x == 0
	}
)

setMethod("isTRUE", signature(x="SpatRaster"),
	function(x) {
		opt <- spatOptions()
		x@ptr <- x@ptr$is_true(opt)
		messages(x, "isTRUE")
	}
)


setMethod("isFALSE", signature(x="SpatRaster"),
	function(x) {
		opt <- spatOptions()
		x@ptr <- x@ptr$is_false(opt)
		messages(x, "isFALSE")
	}
)

setMethod("as.logical", signature(x="SpatRaster"),
	function(x) {
		isTRUE(x)
	}
)


setMethod("is.bool", signature(x="SpatRaster"), 
	function(x) {
		x@ptr$valueType == 3
	}
)
setMethod("is.int", signature(x="SpatRaster"), 
	function(x) {
		x@ptr$valueType == 1
	}
)


setMethod("as.bool", signature(x="SpatRaster"), 
	function(x, filename="", ...) {
		opt <- spatOptions(filename, ...)
		x@ptr <- x@ptr$is_true(opt)
		messages(x, "as.boolean")
	}
)

setMethod("as.int", signature(x="SpatRaster"), 
	function(x, filename="", ...) {
		opt <- spatOptions(filename, ...)
		x@ptr <- x@ptr$math("trunc", opt)
		messages(x, "as.int")
	}
)

setMethod("as.integer", signature(x="SpatRaster"), 
	function(x, filename="", ...) {
		as.int(x, filename, x)
	}
)


setMethod("is.na", signature(x="SpatRaster"),
	function(x) {
		opt <- spatOptions()
		x@ptr <- x@ptr$isnan(opt)
		messages(x, "is.na")
	}
)


setMethod("is.nan", signature(x="SpatRaster"),
	function(x) {
		opt <- spatOptions()
		x@ptr <- x@ptr$isnan(opt)
		messages(x, "is.nan")
	}
)


setMethod("is.finite", signature(x="SpatRaster"),
	function(x) {
		opt <- spatOptions()
		x@ptr <- x@ptr$isfinite(opt)
		messages(x, "is.finite")
	}
)

setMethod("is.infinite", signature(x="SpatRaster"),
	function(x) {
		opt <- spatOptions()
		x@ptr <- x@ptr$isinfinite(opt)
		messages(x, "is.infinite")
	}
)


.summarize <- function(x, ..., fun, na.rm=FALSE, filename="", overwrite=FALSE, wopt=list()) {
	dots <- list(...)
	add <- NULL
	cls <- FALSE
	if (length(dots) > 0) {
		cls <- sapply(dots, function(i) inherits(i, "SpatRaster"))
		if (!all(cls)) {
			add <- dots[!cls]
			if (!is.null(names(add))) {
				error(fun, "additional arguments cannot be names (except for `filename`, `overwrite` and `wopt`)")
			}
			i <- sapply(add, function(x) class(x) %in% c("logical", "integer", "numeric"))
			add <- unlist(add[i], use.names = FALSE)
			if (any(!i)) {
				error(fun, "invalid argument(s)")
			}
		}
	}
	if (any(cls)) {
		x <- sds(c(list(x), dots[cls]))
	} 

	opt <- spatOptions(filename, overwrite, wopt=wopt)
	r <- rast()
	if (is.null(add)) {
		r@ptr <- x@ptr$summary(fun, na.rm, opt)
	} else {
		r@ptr <- x@ptr$summary_numb(fun, add, na.rm, opt)
	}
	messages(r, fun)
	r
}


setMethod("which.max", "SpatRaster",  
	function(x) { 
		opt <- spatOptions()
		x@ptr <- x@ptr$summary("which.max", TRUE, opt)
		messages(x, "which.max")
	}
)

setMethod("which.min", "SpatRaster",  
	function(x) { 
		opt <- spatOptions()
		x@ptr <- x@ptr$summary("which.min", TRUE, opt)
		messages(x, "which.min")
	}
)


setMethod("which.lyr", "SpatRaster",  
	function(x) { 
		opt <- spatOptions()
		x@ptr <- x@ptr$summary("which", TRUE, opt)
		messages(x, "which.lyr")
	}
)



setMethod("Summary", signature(x="SpatRaster"),
	function(x, ..., na.rm=FALSE){
		fun <- as.character(sys.call()[[1L]])
		.summarize(x, ..., fun=fun, na.rm=na.rm)
	}
)

setMethod("Summary", signature(x="SpatVector"),
	function(x, ..., na.rm=FALSE){
		apply(values(x), 2, sys.call()[[1L]], ...)
	}
)




setMethod("Summary", signature(x="SpatExtent"),
	function(x, ..., na.rm=FALSE){
		e <- as.vector(x)
		x <- e[1:2]
		y <- e[3:4]
		fun <- as.character(sys.call()[[1L]])
		if (fun == "range") {
			r <- c(diff(x), diff(y))
			names(r) <- c("x", "y")
			r
		} else {
			c(callGeneric(x), callGeneric(y))
		}
	}
)

setMethod("mean", signature(x="SpatExtent"),
	function(x, ..., trim=NA, na.rm=FALSE){
		if (!is.na(trim)) {	warn("mean", "argument 'trim' is ignored") }
		e <- as.vector(x)
		c(mean(e[1:2]), mean(e[3:4]))
	}
)

setMethod("mean", signature(x="SpatRaster"),
	function(x, ..., trim=NA, na.rm=FALSE){
		if (!is.na(trim)) {	warn("mean", "argument 'trim' is ignored") }
		.summarize(x, ..., fun="mean", na.rm=na.rm)
	}
)

setMethod("mean", signature(x="SpatVector"),
	function(x, ..., trim=NA, na.rm=FALSE){
		if (!is.na(trim)) {	warn("mean", "argument 'trim' is ignored") }
		colMeans(values(x))
	}
)

setMethod("median", signature(x="SpatRaster"),
	function(x, na.rm=FALSE){
		.summarize(x, fun="median", na.rm=na.rm)
	}
)

setMethod("median", signature(x="SpatVector"),
	function(x, na.rm=FALSE){
		apply(values(x), 2, median, na.rm=na.rm)
	}
)

setMethod("Compare", signature(e1="SpatExtent", e2="SpatExtent"),
    function(e1, e2){ 
		oper <- as.vector(.Generic)[1]
		if (!(oper %in% c("==", "!=", ">", "<", ">=", "<="))) {
			error(oper, "is not implemented for SpatExtent")
		}
		return( e1@ptr$compare(e2@ptr, oper, 0.000001) )
	}
)


setMethod("stdev", signature(x="SpatRaster"),
	function(x, ..., pop=TRUE, na.rm=FALSE){
		if (pop) {
			.summarize(x, ..., fun="std", na.rm=na.rm)			
		} else {
			.summarize(x, ..., fun="sd", na.rm=na.rm)
		}
	}
)


setMethod("modal", signature("SpatRaster"), 
	function(x, ..., ties="first", na.rm=FALSE, filename="", overwrite=FALSE, wopt=list()) { 
		opt <- spatOptions(filename, overwrite, wopt=wopt)
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
		messages(x, "modal")
	}
)

