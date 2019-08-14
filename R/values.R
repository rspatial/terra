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
	stopifnot(col > 0)
	stopifnot(row > 0)
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
	if (.hasValues(x)) {
		v <- x@ptr$getValues()
		show_messages(x, "values")
		if (mat) {
			v <- matrix(v, ncol=nlyr(x))
			colnames(v) <- names(x)
		}
	} else {
		v <- NULL
	}
	return(v)	
}
)


setMethod("values<-", signature("SpatRaster", "ANY"), 
	function(x, value) {
	if (is.matrix(value)) { 
		if (nrow(value) == nrow(x)) {
			value <- as.vector(t(value))
		} else {
			value <- as.vector(value)
		} 
	} else if (is.array(value)) { 
		stopifnot(length(dim(value)) == 3)
		value <- as.vector(aperm(value, c(2,1,3)))
	}
	
	if (!(is.numeric(value) || is.integer(value) || is.logical(value))) {
		stop("value must be numeric, integer, or logical")
	}

	lv <- length(value)
	nc <- ncell(x)
	nl <- nlyr(x)
	if (lv == 1) {	
		value <- rep(value, nl * nc)
	} else {
		stopifnot((lv %% nc) == 0)
		if (lv < (nc * nl)) {
			value <- rep(value, length.out=nc*nl)
		}
	}
	y <- rast(x)
	y@ptr$setValues(value)
	y
}
)
	

.hasValues <- function(x) {
	x@ptr$hasValues
}

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
		rbind(rmin, rmax)
	}
)


setMethod("setMinMax", signature(x="SpatRaster"), 
	function(x) {
		if (!.hasMinMax(x)) {
			x@ptr$setRange
		}
	}
)


setMethod("compareGeom", signature(x="SpatRaster", y="SpatRaster"), 
	function(x, y, ..., lyrs=TRUE, crs=TRUE, warncrs=FALSE, ext=TRUE, rowcol=TRUE, res=FALSE) {
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

