# Author: Robert J. Hijmans
# Date :  June 2017
# Version 0.9
# Licence GPL v3


setMethod('dim', signature(x='SpatRaster'), 
	function(x){ return(c(nrow(x), ncol(x), nlyr(x))) }
)

setMethod('nrow', signature(x='SpatRaster'), 
	function(x){ return(x@ptr$nrow)}
)

setMethod('nrow', signature(x='SpatVector'), 
	function(x){ return(x@ptr$nrow())}
)

setMethod('ncol', signature(x='SpatRaster'), 
	function(x){ return(x@ptr$ncol) }
)

setMethod('ncol', signature(x='SpatVector'), 
	function(x){ return(x@ptr$ncol())}
)


setMethod('dim<-', signature(x='SpatRaster'), 
	function(x, value) {
	
		if (length(value) == 1) {
			value <- c(value, ncol(x), nlyr(x))
		} else if (length(value) == 2) {
			value <- c(value, nlyr(x))
		} else if (length(value) > 3) {
			warning('value should have length 1, 2, or 3. Additional values ignored')
			value <- value[1:3]
		}		
		value <- as.integer(pmax(round(value), c(1,1,1)))
		rast(nrow=value[1], ncol=value[2], nlyr=value[3], extent=ext(x), crs=crs(x))
	}
)




setMethod('ncell', signature(x='SpatRaster'), 
	function(x) {
		return(as.numeric(ncol(x)) * nrow(x))
	}
)


setMethod('ncell', signature(x='ANY'), 
	function(x) {
		NROW(x) * NCOL(x)
	}
)



setMethod('length', signature(x='SpatRaster'), 
	function(x) {
		ncell(x) * nlyr(x)
	}
)

if (!isGeneric("nlyr")) {
	setGeneric("nlyr", function(x)
		standardGeneric("nlyr"))
}	

setMethod('nlyr', signature(x='SpatRaster'), 
	function(x){
		return(x@ptr$nlyr() ) 
    }
)


if (!isGeneric("nsrc")) {
	setGeneric("nsrc", function(x)
		standardGeneric("nsrc"))
}	

setMethod('nsrc', signature(x='SpatRaster'), 
	function(x){
		return(x@ptr$nsrc() ) 
    }
)

.nlyrBySource <- function(x) {
	x@ptr$nlyrBySource();
}

'ncol<-' <- function(x, value) {
	dim(x) <- c(nrow(x), value)
	return(x)
}	

'nrow<-' <- function(x, value) {
	dim(x) <- c(value, ncol(x))
	return(x)
}	



setMethod('res', signature(x='SpatRaster'), 
function(x) {
		x@ptr$res
	}
)

setMethod('xres', signature(x='SpatRaster'), 
function(x) {
		res(x)[1]
	}
)

setMethod('yres', signature(x='SpatRaster'), 
function(x) {
		res(x)[2]
	}
)


if (!isGeneric("ext")) {
	setGeneric("ext", function(x, ...)
		standardGeneric("ext"))
}	

if (!isGeneric("ext<-")) {
	setGeneric("ext<-", function(x, value)
		standardGeneric("ext<-"))
}	

