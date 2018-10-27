# Author: Robert J. Hijmans
# Date :  September 2018
# Version 1.0
# Licence GPL v3

# done
# "abs", "sign", "sqrt", "ceiling", "floor", "trunc", "log", "log10", "log2", , 
#"acos", "asin", "atan", "cos", "sin", "tan", 

#to do
# "acosh", "asinh", "atanh", "cosh", "cospi", "sinh", "sinpi", "tanh", "tanpi", 
#"cummax", "cummin", "cumprod", "cumsum",
#"expm1", "log1p"

#"gamma", "lgamma", "digamma", "trigamma")


setMethod("Math", signature(x='SpatRaster'),
    function(x){ 
		oper <- as.vector(.Generic)[1]
		if (oper %in% c("acos", "acosh", "asin", "asinh", "atan", "atanh", "cos", "cosh", "cospi", "sin", "sinh", "sinpi", "tan", "tanh", "tanpi")) {
			x@ptr <- x@ptr$trig(oper, "", FALSE)
		} else {
			x@ptr <- x@ptr$math(oper, "", FALSE)
		}
		.messages(x, oper)
		x
	}	
)

