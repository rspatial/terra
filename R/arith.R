# Author: Robert J. Hijmans
# Date :  October 2018
# Version 1.0
# Licence GPL v3


setMethod("Arith", signature(e1='SpatRaster', e2='SpatRaster'),
    function(e1, e2){ 
		operator <- as.vector(.Generic)[1]
		stopifnot(operator %in% c("+", "-", "*", "/", "%%")) 
		operator <- ifelse(operator == "%%", "%", operator)
		e1@ptr <- e1@ptr$arith_rast(e2@ptr, operator, "", FALSE)
		e1
	}	
)


setMethod("Arith", signature(e1='SpatRaster', e2='numeric'),
    function(e1, e2){ 
		operator <- as.vector(.Generic)[1]
		stopifnot(operator %in% c("+", "-", "*", "/", "%%")) 
		operator <- ifelse(operator == "%%", "%", operator)
		e1@ptr <- e1@ptr$arith_numb(e2, operator, "", FALSE)
		e1
	}	
)



setMethod("Arith", signature(e1='numeric', e2='SpatRaster'),
    function(e1, e2){ 
		operator <- as.vector(.Generic)[1]
		stopifnot(operator %in% c("+", "-", "*", "/", "%%")) 
		operator <- ifelse(operator == "%%", "%", operator)
		e2@ptr <- e2@ptr$arith_rev(e1, operator, "", FALSE)
		e2
	}	
)







setMethod("Compare", signature(e1='SpatRaster', e2='SpatRaster'),
    function(e1, e2){ 
		oper <- as.vector(.Generic)[1]
		e1@ptr <- e1@ptr$arith_rast(e2@ptr, oper, "", FALSE)
		e1
	}	
)


setMethod("Compare", signature(e1='SpatRaster', e2='numeric'),
    function(e1, e2){ 
		oper <- as.vector(.Generic)[1]
		e1@ptr <- e1@ptr$arith_numb(e2, oper, "", FALSE)
		e1
	}
)


setMethod("Compare", signature(e1='numeric', e2='SpatRaster'),
    function(e1, e2){ 
		oper <- as.vector(.Generic)[1]
		e2@ptr <- e2@ptr$arith_numb(e1, oper, "", FALSE)
		e2
	}	
)


