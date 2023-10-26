# Author: Robert J. Hijmans
# Date :  September 2018
# Version 1.0
# License GPL v3

#to do
#"gamma", "lgamma", "digamma", "trigamma")

setMethod("log", signature(x="SpatRaster"),
    function(x, base=exp(1)){
		opt <- spatOptions()
		if (base == exp(1)) {
			x@cpp <- x@cpp$math("log", opt)
		} else if (base == 2) {
			x@cpp <- x@cpp$math("log2", opt)
		} else if (base == 10) {
			x@cpp <- x@cpp$math("log10", opt)
		} else {
			x <- app(x, function(i) log(i, base))
		}
		x
	}
)


#? "gamma", "lgamma", "digamma", "trigamma"
setMethod("Math", signature(x="SpatRaster"),
    function(x){
		oper <- as.vector(.Generic)[1]
		opt <- spatOptions()
		if (substr(oper, 1, 3) == "cum") {
			x@cpp <- x@cpp$cum(substr(oper, 4, 10), FALSE, opt)
		} else if (oper %in% c("acos", "acosh", "asin", "asinh", "atan", "atanh", "cos", "cosh", "cospi", "sin", "sinh", "sinpi", "tan", "tanh", "tanpi")) {
			x@cpp <- x@cpp$trig(oper, opt)
		} else {
			x@cpp <- x@cpp$math(oper, opt)
		}
		messages(x, oper)
	}
)



setMethod("math", signature(x="SpatRaster"),
    function(x, fun, digits=0, filename="", overwrite=FALSE, ...){
		if (!is.character(fun)) {
			error("math", "fun must be a character value")
		}
		fun = fun[1]
		opt <- spatOptions(filename, overwrite, ...)
		if (substr(fun, 1, 3) == "cum") {
			x@cpp <- x@cpp$cum(substr(fun, 4, 10), FALSE, "", FALSE)
		} else if (fun %in% c("acos", "acosh", "asin", "asinh", "atan", "atanh", "cos", "cosh", "cospi", "sin", "sinh", "sinpi", "tan", "tanh", "tanpi")) {
			x@cpp <- x@cpp$trig(fun, opt)
		} else if (fun %in% c("abs", "sign", "sqrt", "ceiling", "floor", "trunc", "log", "log10", "log2")) {		x@cpp <- x@cpp$math(fun, opt)
		} else if (fun %in% c("round", "signif")) {
			x@cpp <- x@cpp$math2(fun, digits, opt)
		} else {
			error("math", "unknown function")
		}
		messages(x, fun)
	}
)



setMethod("Math2", signature(x="SpatRaster"),
    function(x, digits=0){
		opt <- spatOptions()
		oper <- as.vector(.Generic)[1]
		x@cpp <- x@cpp$math2(oper, digits, opt)
		messages(x, oper)
	}
)


setMethod("Math", signature(x="SpatExtent"),
    function(x){
		oper <- as.vector(.Generic)[1]
		if (oper == "floor") {
			x@cpp <- x@cpp$floor()
		} else if (oper == "ceiling") {
			x@cpp <- x@cpp$ceil()
		} else {
			error(oper, "not implemented for SpatExtent")
		}
		if (!is.valid(x)) {
			error(oper, "invalid extent")
		}
		return(x)
	}
)

setMethod("Math2", signature(x="SpatExtent"),
    function(x, digits=0){
		oper <- as.vector(.Generic)[1]
		if (oper == "round") {
			x@cpp <- x@cpp$round(digits)
			if (!is.valid(x)) {
				error(oper, "invalid extent")
			}
			return(x)
		} else {
			error(oper, "not implemented for SpatExtent")
		}
	}
)

setMethod("Math2", signature(x="SpatVector"),
    function(x, digits=4){
		oper <- as.vector(.Generic)[1]
		if (oper == "round") {
			x@cpp <- x@cpp$round(digits)
			return(x)
		} else {
			error(oper, "not implemented for SpatVector")
		}
	}
)

