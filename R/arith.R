# Author: Robert J. Hijmans
# Date :  September 2018
# Version 1.0
# Licence GPL v3


setMethod("Arith", signature(e1='SpatRaster', e2='SpatRaster'),
    function(e1, e2){ 
		oper <- as.vector(.Generic)[1]
		stopifnot(oper %in% c("+", "-", "*", "/", "%%")) 
		oper <- ifelse(oper == "%%", "%", oper)
		e1@ptr <- e1@ptr$arith_rast(e2@ptr, oper, "", FALSE)
		.messages(e1, oper)
		e1
	}	
)


setMethod("Arith", signature(e1='SpatRaster', e2='numeric'),
    function(e1, e2){ 
		oper <- as.vector(.Generic)[1]
		stopifnot(oper %in% c("+", "-", "*", "/", "%%")) 
		oper <- ifelse(oper == "%%", "%", oper)
		e1@ptr <- e1@ptr$arith_numb(e2, oper, "", FALSE)
		.messages(e1, oper)				
		e1
	}	
)


setMethod("Arith", signature(e1='SpatRaster', e2='missing'),
    function(e1, e2){ 
		methods::callGeneric(0, e1)
	}
)

setMethod("Arith", signature(e1='numeric', e2='SpatRaster'),
    function(e1, e2){ 
		oper <- as.vector(.Generic)[1]
		stopifnot(oper %in% c("+", "-", "*", "/", "%%")) 
		oper <- ifelse(oper == "%%", "%", oper)
		e2@ptr <- e2@ptr$arith_rev(e1, oper, "", FALSE)
		.messages(e2, oper)				
		e2
	}	
)




setMethod("Compare", signature(e1='SpatRaster', e2='SpatRaster'),
    function(e1, e2){ 
		oper <- as.vector(.Generic)[1]
		e1@ptr <- e1@ptr$arith_rast(e2@ptr, oper, "", FALSE)
		.messages(e1, oper)
		e1
	}	
)


setMethod("Compare", signature(e1='SpatRaster', e2='numeric'),
    function(e1, e2){ 
		oper <- as.vector(.Generic)[1]
		e1@ptr <- e1@ptr$arith_numb(e2, oper, "", FALSE)
		.messages(e1, oper)
		e1
	}
)


setMethod("Compare", signature(e1='numeric', e2='SpatRaster'),
    function(e1, e2){ 
		oper <- as.vector(.Generic)[1]
		e2@ptr <- e2@ptr$arith_rev(e1, oper, "", FALSE)
		.messages(e2, oper)
		e2
	}	
)


setMethod("Logic", signature(e1='SpatRaster', e2='SpatRaster'),
    function(e1, e2){ 
		oper <- as.vector(.Generic)[1]
		e1@ptr <- e1@ptr$logic(e2@ptr, oper, "", FALSE)
		.messages(e1, oper)
		e1
	}	
)

setMethod("Logic", signature(e1='SpatRaster', e2='numeric'),
    function(e1, e2){ 
		oper <- as.vector(.Generic)[1]
		e1@ptr <- e1@ptr$logic_numb(e2, oper, "", FALSE)
		.messages(e1, oper)
		e1
	}
)


setMethod("Logic", signature(e1='numeric', e2='SpatRaster'),
    function(e1, e2){ 
		oper <- as.vector(.Generic)[1]
		e2@ptr <- e2@ptr$logic_numb(e1, oper, "", FALSE)
		.messages(e2, oper)
		e2
	}	
)
