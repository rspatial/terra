# Author: Robert J. Hijmans
# Date :  June 2017
# Version 0.9
# License GPL v3


setMethod("dim", signature(x="SpatRaster"), 
	function(x){ return(c(nrow(x), ncol(x), nlyr(x))) }
)

setMethod("dim", signature(x="SpatRasterDataset"), 
	function(x) {
		dim(x[1])[1:2]
	}
)

setMethod("nrow", signature(x="SpatRaster"), 
	function(x){ return(x@ptr$nrow())}
)

setMethod("nrow", signature(x="SpatRasterDataset"), 
	function(x){ return(x[1]@ptr$nrow())}
)

setMethod("nrow", signature(x="SpatVector"), 
	function(x){ return(x@ptr$nrow())}
)

setMethod("ncol", signature(x="SpatRaster"), 
	function(x){ return(x@ptr$ncol()) }
)

setMethod("ncol", signature(x="SpatRasterDataset"), 
	function(x){ return(x[1]@ptr$ncol())}
)

setMethod("ncol", signature(x="SpatVector"), 
	function(x){ return(x@ptr$ncol())}
)


setMethod("dim<-", signature(x="SpatRaster"), 
	function(x, value) {
	
		if (length(value) == 1) {
			value <- c(value, ncol(x), nlyr(x))
		} else if (length(value) == 2) {
			value <- c(value, nlyr(x))
		} else if (length(value) > 3) {
			warning("value should have length 1, 2, or 3. Additional values ignored")
			value <- value[1:3]
		}		
		value <- as.integer(pmax(round(value), c(1,1,1)))
		rast(nrow=value[1], ncol=value[2], nlyr=value[3], extent=ext(x), crs=crs(x))
	}
)



setMethod("ncell", signature(x="SpatRaster"), 
	function(x) {
		return(as.numeric(ncol(x)) * nrow(x))
	}
)

setMethod("ncell", signature(x="SpatRasterDataset"), 
	function(x) {
		ncell(x[1])
	}
)


setMethod("ncell", signature(x="ANY"), 
	function(x) {
		NROW(x) * NCOL(x)
	}
)


setMethod("size", signature(x="SpatRaster"), 
	function(x) {
		x@ptr$size()
	}
)

setMethod("size", signature(x="SpatRasterDataset"), 
	function(x){
		nc <- ncell(x)
		sapply(1:length(x), function(i) nlyr(x[i]) * nc)
    }
)

setMethod("nlyr", signature(x="SpatRaster"), 
	function(x){
		return(x@ptr$nlyr() ) 
    }
)

setMethod("nlyr", signature(x="SpatRasterDataset"), 
	function(x){
		sapply(1:length(x), function(i) nlyr(x[i]))
    }
)


setMethod("nsrc", signature(x="SpatRaster"), 
	function(x){
		return(x@ptr$nsrc() ) 
    }
)

.nlyrBySource <- function(x) {
	x@ptr$nlyrBySource();
}


setMethod("ncol<-", signature("SpatRaster", "numeric"), 
	function(x, ..., value) {
		dim(x) <- c(nrow(x), value)
		return(x)
	}
)

setMethod("nrow<-", signature("SpatRaster", "numeric"), 
	function(x, ..., value) {
		dim(x) <- c(value, ncol(x))
		return(x)
	}
)

setMethod("nlyr<-", signature("SpatRaster", "numeric"), 
	function(x, ..., value) {
		dim(x) <- c(nrow(x), ncol(x), value)
		return(x)
	}
)


setMethod("res", signature(x="SpatRaster"), 
function(x) {
		x@ptr$res
	}
)

setMethod("res", signature(x="SpatRasterDataset"), 
function(x) {
		x[1]@ptr$res
	}
)

setMethod("res<-", signature(x="SpatRaster"), 
	function(x, value) {
		if (length(value) == 1) {
			value <- c(value, value)
		} else if (length(value) > 2) {
			warning("value should have length 1 or 2. Additional values ignored")
		}		
		x@ptr <- x@ptr$set_resolution(value[1], value[2])
		show_messages(x, "resolution")
	}
)


setMethod("xres", signature(x="SpatRaster"), 
function(x) {
		res(x)[1]
	}
)

setMethod("yres", signature(x="SpatRaster"), 
function(x) {
		res(x)[2]
	}
)


