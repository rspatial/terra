# Author: Robert J. Hijmans
# Date :  June 2008
# Version 0.9
# License GPL v3

setMethod("hasValues", signature(x="SpatRaster"), 
	function(x, ...) {
		x@ptr$hasValues
	}
)

setMethod("readValues", signature(x="SpatRaster"), 
function(x, row=1, nrows=nrow(x), col=1, ncols=ncol(x), mat=FALSE, dataframe=FALSE, ...) {
	stopifnot(row > 0)
	stopifnot(col > 0)
	v <- x@ptr$readValues(row-1, nrows, col-1, ncols)
	if (dataframe || mat) {
		v <- matrix(v, ncol = nlyr(x))
		colnames(v) <- names(x)	
	}
	if (dataframe) {
		v <- data.frame(v)	
		ff <- is.factor(x)
		if (any(ff)) {
			ff <- which(ff)
			levs <- levels(x)
			for (f in ff) {
				lev <- levs[[f]]
				v[[f]] = factor(v[[f]], levels=lev$levels)
				levels(v[[f]]) = lev$labels
			}
		}
	}
	v
}
)


setMethod("values", signature(x="SpatRaster"), 
function(x, mat=TRUE, ...) {
	if (hasValues(x)) {
		v <- x@ptr$getValues(-1)
		show_messages(x, "values")
	} else {
		v <- matrix(NA, nrow=ncell(x), ncol=nlyr(x))
	}
	if (mat) {
		v <- matrix(v, ncol=nlyr(x))
		colnames(v) <- names(x)
	}
	return(v)	
}
)

setMethod("values<-", signature("SpatRaster", "ANY"), 
	function(x, value) {
		setValues(x, value)
	}
)

setMethod("setValues", signature("SpatRaster", "ANY"), 
	function(x, values, ...) {

		if (is.matrix(values)) { 
			if (nrow(values) == nrow(x)) {
				values <- as.vector(t(values))
			} else {
				values <- as.vector(values)
			} 
		} else if (is.array(values)) { 
			stopifnot(length(dim(values)) == 3)
			values <- as.vector(aperm(values, c(2,1,3)))
		}
		
		if (!(is.numeric(values) || is.integer(values) || is.logical(values))) {
			stop("values must be numeric, integer, or logical")
		}

		lv <- length(values)
		nc <- ncell(x)
		nl <- nlyr(x)
		if (lv == 1) {	
			values <- rep(values, nl * nc)
		} else {
			if (!((lv %% nc) == 0)) {
				warning("the length of the values does not match the size of the SpatRaster")
			}
			if (lv > (nc * nl)) {
				values <- values[1:(nc*nl)]
			} else if (lv < (nc * nl)) {
				values <- rep(values, length.out=nc*nl)
			}
		}
		y <- rast(x)
		y@ptr$setValues(values)
		y
	}
)
	

#.hasValues <- function(x) {
#	x@ptr$hasValues
#}

.inMemory <- function(x) {
	x@ptr$inMemory
}

.filenames <- function(x) {
	x@ptr$filenames
}

.hasMinMax <- function(x) {
	x@ptr$hasRange
}



setMethod("sources", signature(x="SpatRaster"), 
	function(x, ...) {
		src <- x@ptr$filenames
		src[src == ""] <= "memory"
		data.frame(source=src, nlyr=x@ptr$nlyrBySource(), stringsAsFactors=FALSE)
	}
)


setMethod("minmax", signature(x="SpatRaster"), 
	function(x) {
		rmin <- x@ptr$range_min
		rmax <- x@ptr$range_max
		r <- rbind(rmin, rmax)
		colnames(r) <- names(x)
		r
	}
)


setMethod("setMinMax", signature(x="SpatRaster"), 
	function(x) {
		if (any(!.hasMinMax(x))) {
			x@ptr$setRange()
			x <- show_messages(x)
		}
	}
)


setMethod("compareGeom", signature(x="SpatRaster", y="SpatRaster"), 
	function(x, y, ..., lyrs=FALSE, crs=TRUE, warncrs=FALSE, ext=TRUE, rowcol=TRUE, res=FALSE) {
		dots <- list(...)
		bool <- x@ptr$compare_geom(y@ptr, lyrs, crs, warncrs, ext, rowcol, res)
		show_messages(x, "compareGeom")
		if (length(dots)>1) {
			for (i in 1:length(dots)) {
				bool <- x@ptr$compare_geom(dots[[i]]@ptr, lyrs, crs, warncrs, ext, rowcol, res)
				show_messages(x, "compareGeom")
			}
		}
		bool
	}
)


setMethod("values", signature("SpatVector"), 
	function(x, ...) {
		as.data.frame(x, ...)
	}
)


setMethod("values<-", signature("SpatVector", "data.frame"), 
	function(x, value) {
		stopifnot(nrow(x) == nrow(value))
		x <- x[,0]
		types <- sapply(value, class)
		nms <- colnames(value)
		for (i in 1:ncol(value)) {
			if (types[i] == "numeric") {
				x@ptr$add_column_double(value[[i]], nms[i])
			} else if (types[i] == "integer") {
				x@ptr$add_column_long(value[[i]], nms[i])
			} else if (types[i] == "character") {
				x@ptr$add_column_string(value[[i]], nms[i])
			} else {
				att <- as.character(value[[i]])
				x@ptr$add_column_string(att, nms[i])
			}
		}
		x
	}
)

