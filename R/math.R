# Author: Robert J. Hijmans
# Date :  September 2018
# Version 1.0
# Licence GPL v3

#to do
#"gamma", "lgamma", "digamma", "trigamma")


setMethod("Math", signature(x='SpatRaster'),
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

