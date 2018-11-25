# Author: Robert J. Hijmans
# Date :  September 2018
# Version 1.0
# License GPL v3


setMethod("Arith", signature(e1='SpatRaster', e2='SpatRaster'),
    function(e1, e2){ 
		oper <- as.vector(.Generic)[1]
		stopifnot(oper %in% c("+", "-", "*", "/", "%%")) 
		oper <- ifelse(oper == "%%", "%", oper)
		e1@ptr <- e1@ptr$arith_rast(e2@ptr, oper, .terra_environment$options@ptr)
		show_messages(e1, oper)
	}	
)


setMethod("Arith", signature(e1='SpatRaster', e2='numeric'),
    function(e1, e2){ 
		oper <- as.vector(.Generic)[1]
		stopifnot(oper %in% c("+", "-", "*", "/", "%%")) 
		oper <- ifelse(oper == "%%", "%", oper)
		e1@ptr <- e1@ptr$arith_numb(e2, oper, .terra_environment$options@ptr)
		show_messages(e1, oper)				
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
		e2@ptr <- e2@ptr$arith_rev(e1, oper, .terra_environment$options@ptr)
		show_messages(e2, oper)				
	}	
)




setMethod("Compare", signature(e1='SpatRaster', e2='SpatRaster'),
    function(e1, e2){ 
		oper <- as.vector(.Generic)[1]
		e1@ptr <- e1@ptr$arith_rast(e2@ptr, oper, .terra_environment$options@ptr)
		show_messages(e1, oper)
	}	
)


setMethod("Compare", signature(e1='SpatRaster', e2='numeric'),
    function(e1, e2){ 
		oper <- as.vector(.Generic)[1]
		e1@ptr <- e1@ptr$arith_numb(e2, oper, .terra_environment$options@ptr)
		show_messages(e1, oper)
	}
)


setMethod("Compare", signature(e1='numeric', e2='SpatRaster'),
    function(e1, e2){ 
		oper <- as.vector(.Generic)[1]
		e2@ptr <- e2@ptr$arith_rev(e1, oper, .terra_environment$options@ptr)
		show_messages(e2, oper)
	}	
)


setMethod("Logic", signature(e1='SpatRaster', e2='SpatRaster'),
    function(e1, e2){ 
		oper <- as.vector(.Generic)[1]
		e1@ptr <- e1@ptr$logic(e2@ptr, oper, .terra_environment$options@ptr)
		show_messages(e1, oper)
	}	
)

setMethod("Logic", signature(e1='SpatRaster', e2='numeric'),
    function(e1, e2){ 
		oper <- as.vector(.Generic)[1]
		e1@ptr <- e1@ptr$logic_numb(e2, oper, .terra_environment$options@ptr)
		show_messages(e1, oper)
	}
)


setMethod("Logic", signature(e1='numeric', e2='SpatRaster'),
    function(e1, e2){ 
		oper <- as.vector(.Generic)[1]
		e2@ptr <- e2@ptr$logic_numb(e1, oper, .terra_environment$options@ptr)
		show_messages(e2, oper)
	}	
)
