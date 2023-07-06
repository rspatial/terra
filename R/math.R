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
			x@pnt <- x@pnt$math("log", opt)
		} else if (base == 2) {
			x@pnt <- x@pnt$math("log2", opt)
		} else if (base == 10) {
			x@pnt <- x@pnt$math("log10", opt)
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
			x@pnt <- x@pnt$cum(substr(oper, 4, 10), FALSE, opt)
		} else if (oper %in% c("acos", "acosh", "asin", "asinh", "atan", "atanh", "cos", "cosh", "cospi", "sin", "sinh", "sinpi", "tan", "tanh", "tanpi")) {
			x@pnt <- x@pnt$trig(oper, opt)
		} else {
			x@pnt <- x@pnt$math(oper, opt)
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
			x@pnt <- x@pnt$cum(substr(fun, 4, 10), FALSE, "", FALSE)
		} else if (fun %in% c("acos", "acosh", "asin", "asinh", "atan", "atanh", "cos", "cosh", "cospi", "sin", "sinh", "sinpi", "tan", "tanh", "tanpi")) {
			x@pnt <- x@pnt$trig(fun, opt)
		} else if (fun %in% c("abs", "sign", "sqrt", "ceiling", "floor", "trunc", "log", "log10", "log2")) {		x@pnt <- x@pnt$math(fun, opt)
		} else if (fun %in% c("round", "signif")) {
			x@pnt <- x@pnt$math2(fun, digits, opt)
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
		x@pnt <- x@pnt$math2(oper, digits, opt)
		messages(x, oper)
	}
)


setMethod("Math", signature(x="SpatExtent"),
    function(x){
		oper <- as.vector(.Generic)[1]
		if (oper == "floor") {
			x@pnt <- x@pnt$floor()
		} else if (oper == "ceiling") {
			x@pnt <- x@pnt$ceil()
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
			x@pnt <- x@pnt$round(digits)
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
			x@pnt <- x@pnt$round(digits)
			return(x)
		} else {
			error(oper, "not implemented for SpatVector")
		}
	}
)

