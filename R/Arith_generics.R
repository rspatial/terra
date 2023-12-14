# Author: Robert J. Hijmans
# Date :  September 2018
# Version 1.0
# License GPL v3


setMethod("Arith", signature(e1="SpatExtent", e2="numeric"),
    function(e1, e2){
		oper <- as.vector(.Generic)[1]
		if (oper == "%%") {
			e1@cpp <- e1@cpp$align(e2[1], "")
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
		if (!is.valid(e1)) {
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
		if (!is.valid(e1)) {
			error(oper, "this would create an invalid extent")
		}
		e1
	}
)


setMethod("Arith", signature(e1="SpatExtent", e2="SpatExtent"),
    function(e1, e2){
		oper <- as.vector(.Generic)[1]
		if (oper == "+") {
			e1@cpp <- e1@cpp$deepcopy()
			e1@cpp$union(e2@cpp)
		} else if (oper == "*") {
			e1@cpp <- e1@cpp$intersect(e2@cpp)
		} else if (oper == "/") {
			d <- c(diff(e1[1:2]) / diff(e2[1:2]), diff(e1[3:4]) / diff(e2[3:4]))
			names(d) <- c("x", "y")
			return(d)
		} else {
			error(oper, "only +, *, and / are supported for SpatExtent")
		}
		if (!is.valid(e1)) {
			error(oper, "this would create an invalid extent")
		}
		e1
	}
)



setMethod("Arith", signature(e1="SpatVector", e2="SpatVector"),
    function(e1, e2){
		oper <- as.vector(.Generic)[1]
		if (geomtype(e1) != geomtype(e2)) {
			error(oper, "geometry types do not match")
		}
		if (oper == "+") {
			e1 <- union(e1, e2)
		} else if (oper == "*") {
			e1 <- intersect(e1, e2)
		} else if (oper == "-") {
			e1 <- erase(e1, e2)
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
		e1@cpp <- e1@cpp$arith_rast(e2@cpp, oper, FALSE, opt)
		messages(e1, oper)
	}
)


setMethod("Arith", signature(e1="SpatRaster", e2="numeric"),
    function(e1, e2){
		oper <- as.vector(.Generic)[1]
		opt <- spatOptions()
		e1@cpp <- e1@cpp$arith_numb(e2, oper, FALSE, FALSE, opt)
		messages(e1, oper)
	}
)

setMethod("Arith", signature(e1="SpatRaster", e2="logical"),
    function(e1, e2){
		methods::callGeneric(e1, as.integer(e2))
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
		opt <- spatOptions()
		e2@cpp <- e2@cpp$arith_numb(e1, oper, TRUE, FALSE, opt)
		messages(e2, oper)
	}
)

setMethod("Arith", signature(e1="logical", e2="SpatRaster"),
    function(e1, e2){
		methods::callGeneric(as.integer(e1), e2)
	}
)

setMethod("Arith", signature(e1="SpatRaster", e2="matrix"),
    function(e1, e2){
		oper <- as.vector(.Generic)[1]
		opt <- spatOptions()
		e1@cpp <- e1@cpp$arith_m(as.vector(e2), oper, dim(e2)[1:2], FALSE, opt)
		messages(e1, oper)
	}
)

setMethod("Arith", signature(e1="matrix", e2="SpatRaster"),
    function(e1, e2){
		oper <- as.vector(.Generic)[1]
		opt <- spatOptions()
		e1@cpp <- e1@cpp$arith_m(as.vector(e2), oper, dim(e2)[1:2], TRUE, opt)
		messages(e1, oper)
	}
)



setMethod("Compare", signature(e1="SpatRaster", e2="SpatRaster"),
    function(e1, e2){
		opt <- spatOptions()
		oper <- as.vector(.Generic)[1]
		e1@cpp <- e1@cpp$arith_rast(e2@cpp, oper, FALSE, opt)
		messages(e1, oper)
	}
)


setMethod("Compare", signature(e1="SpatRaster", e2="numeric"),
    function(e1, e2){
		oper <- as.vector(.Generic)[1]
		opt <- spatOptions()
		e1@cpp <- e1@cpp$arith_numb(e2, oper, FALSE, FALSE, opt)
		messages(e1, oper)
	}
)


setMethod("Compare", signature(e1="numeric", e2="SpatRaster"),
    function(e1, e2){
		oper <- as.vector(.Generic)[1]
		opt <- spatOptions()
		e2@cpp <- e2@cpp$arith_numb(e1, oper, TRUE, FALSE, opt)
		messages(e2, oper)
	}
)

setMethod("Compare", signature(e1="SpatRaster", e2="matrix"),
    function(e1, e2){
		oper <- as.vector(.Generic)[1]
		opt <- spatOptions()
		e1@cpp <- e1@cpp$arith_m(as.vector(e2), oper, dim(e2)[1:2], FALSE, opt)
		messages(e1, oper)
	}
)

setMethod("Compare", signature(e1="matrix", e2="SpatRaster"),
    function(e1, e2){
		oper <- as.vector(.Generic)[1]
		opt <- spatOptions()
		e1@cpp <- e1@cpp$arith_m(as.vector(e2), oper, dim(e2)[1:2], TRUE, opt)
		messages(e1, oper)
	}
)


getFactTable <- function(x, table, sender="%in%") {
	if (!is.factor(x)) {
		error(sender, "Can only match character values if x is categorical")
	}
	if (nlyr(x) != 1) {
		error(sender, "matching with character values is only supported for single layer SpatRaster")
	}
	d <- levels(x)[[1]]
	m <- stats::na.omit(match(table, d[,2]))
	if (length(m) == 0) {
		return(as.logical(x*0))
	}
	d[m,1]
}

setMethod("Compare", signature(e1="SpatRaster", e2="character"),
    function(e1, e2){
		oper <- as.vector(.Generic)[1]
		e2 <- getCatIDs(e1, e2, "==")
		if (oper != "==") {
			error(oper, "only '==' is supported with categorical comparisons")
		}
		if (length(e2) == 0) {
			return(as.logical(e1*0))
		}
		if (length(e2) != 1) {
			error(oper, "comparisons only supported for single values (see %in% and match)")
		}
		opt <- spatOptions()
		e1@cpp <- e1@cpp$arith_numb(e2, oper, TRUE, FALSE, opt)
		messages(e1, oper)
	}
)



setMethod("Logic", signature(e1="SpatRaster", e2="SpatRaster"),
    function(e1, e2){
		oper <- as.vector(.Generic)[1]
		opt <- spatOptions()
		e1@cpp <- e1@cpp$logic_rast(e2@cpp, oper, opt)
		messages(e1, oper)
	}
)

setMethod("Logic", signature(e1="SpatRaster", e2="numeric"),
    function(e1, e2){
		opt <- spatOptions()
		oper <- as.vector(.Generic)[1]
		e2 <- as.logical(e2)
		e1@cpp <- e1@cpp$logic_numb(e2, oper, opt)
		messages(e1, oper)
	}
)


setMethod("Logic", signature(e1="numeric", e2="SpatRaster"),
    function(e1, e2){
		opt <- spatOptions()
		oper <- as.vector(.Generic)[1]
		e1 <- as.logical(e1)
		e2@cpp <- e2@cpp$logic_numb(e1, oper, opt)
		messages(e2, oper)
	}
)

setMethod("Logic", signature(e1="SpatRaster", e2="logical"),
    function(e1, e2){
		opt <- spatOptions()
		oper <- as.vector(.Generic)[1]
		e1@cpp <- e1@cpp$logic_numb(e2, oper, opt)
		messages(e1, oper)
	}
)


setMethod("Logic", signature(e1="logical", e2="SpatRaster"),
    function(e1, e2){
		opt <- spatOptions()
		oper <- as.vector(.Generic)[1]
		e2@cpp <- e2@cpp$logic_numb(e1, oper, opt)
		messages(e2, oper)
	}
)


setMethod("!", signature(x="SpatRaster"),
	function(x) {
		x == 0
	}
)

setMethod("not.na", signature(x="SpatRaster"),
	function(x, falseNA=FALSE, filename="", ...) {
		opt <- spatOptions(filename=filename, ...)
		x@cpp <- x@cpp$not_na(falseNA, opt)
		messages(x, "not.na")
	}
)


setMethod("isTRUE", signature(x="SpatRaster"),
	function(x) {
		opt <- spatOptions()
		x@cpp <- x@cpp$is_true(FALSE, opt)
		messages(x, "isTRUE")
	}
)


setMethod("isFALSE", signature(x="SpatRaster"),
	function(x) {
		opt <- spatOptions()
		x@cpp <- x@cpp$is_false(FALSE, opt)
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
		x@cpp$valueType(FALSE) == 3
	}
)
setMethod("is.int", signature(x="SpatRaster"),
	function(x) {
		x@cpp$valueType(FALSE) == 1
	}
)


setMethod("as.bool", signature(x="SpatRaster"),
	function(x, filename="", ...) {
		opt <- spatOptions(filename, ...)
		x@cpp <- x@cpp$is_true(FALSE, opt)
		messages(x, "as.boolean")
	}
)

setMethod("as.int", signature(x="SpatRaster"),
	function(x, filename="", ...) {
		opt <- spatOptions(filename, ...)
		x@cpp <- x@cpp$math("trunc", opt)
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
		x@cpp <- x@cpp$isnan(FALSE, opt)
		messages(x, "is.na")
	}
)

setMethod("countNA", signature(x="SpatRaster"),
	function(x, n=0) {
		opt <- spatOptions()
		n <- round(n)
		if (n == 1) {
			x@cpp <- x@cpp$anynan(FALSE, opt)
		} else {
			x@cpp <- x@cpp$countnan(n, opt)
		}
		messages(x, "countNA")
	}
)


setMethod("anyNA", signature(x="SpatRaster"),
	function(x) {
		opt <- spatOptions()
		x@cpp <- x@cpp$anynan(FALSE, opt)
		messages(x, "anyNA")
	}
)

setMethod("noNA", signature(x="SpatRaster"),
	function(x, falseNA=FALSE) {
		opt <- spatOptions()
		x@cpp <- x@cpp$nonan(falseNA, opt)
		messages(x, "noNA")
	}
)


setMethod("allNA", signature(x="SpatRaster"),
	function(x, falseNA=FALSE) {
		opt <- spatOptions()
		x@cpp <- x@cpp$allnan(falseNA, opt)
		messages(x, "allNA")
	}
)

setMethod("is.nan", signature(x="SpatRaster"),
	function(x) {
		opt <- spatOptions()
		x@cpp <- x@cpp$isnan(FALSE, opt)
		messages(x, "is.nan")
	}
)


setMethod("is.finite", signature(x="SpatRaster"),
	function(x) {
		opt <- spatOptions()
		x@cpp <- x@cpp$isfinite(FALSE, opt)
		messages(x, "is.finite")
	}
)

setMethod("is.infinite", signature(x="SpatRaster"),
	function(x) {
		opt <- spatOptions()
		x@cpp <- x@cpp$isinfinite(FALSE, opt)
		messages(x, "is.infinite")
	}
)


.summarize <- function(x, ..., fun, na.rm=FALSE, filename="", overwrite=FALSE, wopt=list(), par=FALSE) {

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
	if (any(cls) | par) {
		x <- sds(c(list(x), dots[cls]))
	}

	opt <- spatOptions(filename, overwrite, wopt=wopt)
	r <- rast()
	if (is.null(add)) {
		r@cpp <- x@cpp$summary(fun, na.rm, opt)
	} else {
		r@cpp <- x@cpp$summary_numb(fun, add, na.rm, opt)
	}
	messages(r, fun)
	r
}


setMethod("which.max", "SpatRaster",
	function(x) {
		opt <- spatOptions()
		x@cpp <- x@cpp$summary("which.max", TRUE, opt)
		messages(x, "which.max")
	}
)

setMethod("which.min", "SpatRaster",
	function(x) {
		opt <- spatOptions()
		x@cpp <- x@cpp$summary("which.min", TRUE, opt)
		messages(x, "which.min")
	}
)


setMethod("which.lyr", "SpatRaster",
	function(x) {
		opt <- spatOptions()
		x@cpp <- x@cpp$summary("which", TRUE, opt)
		messages(x, "which.lyr")
	}
)

wherefun <- function(out, list, values) {
	if (list) {
		if (values) {
			lapply(out, function(i) {
				m <- matrix(i, ncol=2)
				m[,1] <- m[,1] + 1
				colnames(m) <- c("cell", "value")
				m
			})
		} else {
			lapply(out, function(i) {i + 1})
		}
	} else {
		if (values) {
			out <- lapply(1:length(out), function(i) {
				m <- matrix(out[[i]], ncol=2)
				m[,1] <- m[,1] + 1
				cbind(i, m)
			})
			out <- do.call(rbind, out)
			colnames(out) <- c("layer", "cell", "value")
		} else {
			out <- lapply(1:length(out), function(i) {cbind(i, out[[i]] + 1)})
			out <- do.call(rbind, out)
			colnames(out) <- c("layer", "cell")
		}
		out
	}
}


setMethod("where.max", "SpatRaster",
	function(x, values=TRUE, list=FALSE) {
		opt <- spatOptions()
		out <- x@cpp$where("max", values, opt)
		x <- messages(x, "where.max")
		wherefun(out, list, values)
	}
)

setMethod("where.min", "SpatRaster",
	function(x, values=TRUE, list=FALSE) {
		opt <- spatOptions()
		out <- x@cpp$where("min", values, opt)
		x <- messages(x, "where.min")
		wherefun(out, list, values)
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
		if (!is.null(list(...))) {	warn("mean", "additional arguments are ignored") }
		colMeans(values(x))
	}
)

setMethod("median", signature(x="SpatRaster"),
	function(x, na.rm=FALSE, ...){
		if (!is.logical(na.rm)) {
			error("median", "na.rm (the second argument) must be a logical value")
		}
		.summarize(x, ..., fun="median", na.rm=na.rm)
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
		return( e1@cpp$compare(e2@cpp, oper, 0.000001) )
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
		x@cpp <- x@cpp$modal(add, ties[1], na.rm[1], opt)
		messages(x, "modal")
	}
)


setMethod("compare", signature(x="SpatRaster"),
    function(x, y, oper, falseNA=FALSE, filename="", overwrite=FALSE, ...){
		if (!is.character(oper)) {
			error("compare", "oper must be a character value")
		}
		oper = oper[1]
		ops <- c("==", "!=", ">", "<", ">=", "<=")
		if (!(oper %in% ops)) {
			error("compare", "oper must be a one of", ops)		
		}

		opt <- spatOptions(filename, overwrite, ...)
		if (inherits(y, "SpatRaster")) {
			x@cpp <- x@cpp$arith_rast(y@cpp, oper, falseNA[1], opt)
		} else {
			x@cpp <- x@cpp$arith_numb(y, oper, FALSE, falseNA[1], opt)
		}
		messages(x, oper)
	}
)


setMethod("logic", signature(x="SpatRaster"),
    function(x, oper, falseNA=FALSE, filename="", overwrite=FALSE, ...){
		if (!is.character(oper)) {
			error("logic", "oper must be a character value")
		}
		oper = oper[1]
		ops <- c("!", "is.na", "allNA", "noNA", "is.infinite", "is.finite", "iSTRUE", "isFALSE")
		if (!(oper %in% ops)) {
			error("compare", "oper must be a one of", ops)		
		}
		opt <- spatOptions(filename, overwrite, ...)
		falseNA <- as.logical(falseNA[1])
		
		if (oper == "is.infinite") {
			x@cpp <- x@cpp$isinfinite(falseNA, opt)
		} else if (oper == "is.finite") {
			x@cpp <- x@cpp$isfinite(falseNA, opt)
		} else if (oper == "is.na") {
			x@cpp <- x@cpp$isnan(falseNA, opt)
		} else if (oper == "isTRUE") {
			x@cpp <- x@cpp$is_true(falseNA, opt)
		} else if (oper == "isFALSE") {
			x@cpp <- x@cpp$is_false(falseNA, opt)
		} else if (oper == "allNA") {
			x@cpp <- x@cpp$allnan(falseNA, opt)
		} else if (oper == "noNA") {
			x@cpp <- x@cpp$nonan(falseNA, opt)
		} else if (oper == "anyNA") {
			x@cpp <- x@cpp$anynan(falseNA, opt)
		} else if (oper == "anyNA") {
			x@cpp <- x@cpp$anynan(falseNA, opt)
		} else if (oper == "!") {
			x@cpp <- x@cpp$arith_numb(0, "==", FALSE, falseNA[1], opt)
		} else {
			error("logic", "??")
		}
		messages(x, "logic")
	}
)

