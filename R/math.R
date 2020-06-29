# Author: Robert J. Hijmans
# Date :  September 2018
# Version 1.0
# License GPL v3

#to do
#"gamma", "lgamma", "digamma", "trigamma")


setMethod("Math", signature(x="SpatRaster"),
    function(x){ 
		oper <- as.vector(.Generic)[1]
		if (substr(oper, 1, 3) == "cum") {
			x@ptr <- x@ptr$cum(substr(oper, 4, 10), FALSE, "", FALSE)
		} else if (oper %in% c("acos", "acosh", "asin", "asinh", "atan", "atanh", "cos", "cosh", "cospi", "sin", "sinh", "sinpi", "tan", "tanh", "tanpi")) {
			x@ptr <- x@ptr$trig(oper, .terra_environment$options@ptr)
		} else {
			x@ptr <- x@ptr$math(oper, .terra_environment$options@ptr)
		}
		show_messages(x, oper)
	}	
)


setMethod("Math2", signature(x="SpatRaster"),
    function(x, digits=0){ 
		oper <- as.vector(.Generic)[1]
		x@ptr <- x@ptr$math2(oper, digits, .terra_environment$options@ptr)
		show_messages(x, oper)
	}	
)


setMethod("Math", signature(x="SpatExtent"),
    function(x){ 
		oper <- as.vector(.Generic)[1]
		if (oper == "floor") {
			x@ptr <- x@ptr$floor()
		} else if (oper == "ceiling") {
			x@ptr <- x@ptr$ceil()
		} else {
			stop("not implemented for SpatExtent")
		}
		if (!x@ptr$valid) {
			stop("invalid extent")
		}
		return(x)		
	}	
)

setMethod("Math2", signature(x="SpatExtent"),
    function(x, digits=0){ 
		oper <- as.vector(.Generic)[1]
		if (oper == "round") {
			x@ptr <- x@ptr$round(x@ptr, digits)
			if (!x@ptr$valid) {
				stop("invalid extent")
			}
			return(x)
		} else {
			stop("not implemented for SpatExtent")
		}
	}	
)


