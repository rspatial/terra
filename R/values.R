# Author: Robert J. Hijmans
# Date :  June 2008
# Version 0.9
# License GPL v3

setMethod("hasValues", signature(x="SpatRaster"), 
	function(x) {
		x@ptr$hasValues
	}
)

setMethod("readValues", signature(x="SpatRaster"), 
function(x, row=1, nrows=nrow(x), col=1, ncols=ncol(x), mat=FALSE, dataframe=FALSE) {
	stopifnot(row > 0 && nrows > 0)
	stopifnot(col > 0 && ncols > 0)
	
	v <- x@ptr$readValues(row-1, nrows, col-1, ncols)
	messages(x, "readValues")
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
				fct <- levs[[f]]
				v[[f]] = factor(v[[f]], levels=(1:length(fct))-1)
				levels(v[[f]]) = fct
			}
		}
	}
	v
}
)


setMethod("values", signature(x="SpatRaster"), 
function(x, mat=TRUE, dataframe=FALSE, row=1, nrows=nrow(x), col=1, ncols=ncol(x)) {
	readStart(x)
	on.exit(readStop(x))
	v <- readValues(x, row, nrows, col, ncols, mat=mat, dataframe=dataframe)
	messages(x)
	return(v)
}
)

setMethod("values<-", signature("SpatRaster", "ANY"), 
	function(x, value) {
		setValues(x, value)
	}
)

setMethod("focalValues", signature("SpatRaster"), 
	function(x, w=3, row=1, nrows=nrow(x), fill=NA) {
		if (is.matrix(w)) {
			#m <- as.vector(t(w))
			w <- dim(w)
		} else {
			w <- rep_len(w, 2)
		}
		readStart(x)
		on.exit(readStop(x))
		m <- matrix(x@ptr$focalValues(w, fill, row-1, nrows), ncol=prod(w), byrow=TRUE)
		messages(x)
		m
	}
)


setMethod("setValues", signature("SpatRaster", "ANY"), 
	function(x, values) {

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
		make_factor <- FALSE
		if (is.character(values)) {
			values <- as.factor(values)
			levs <- levels(values)
			values <- as.integer(values) - 1
			if (max(values, na.rm=TRUE) <= 255) {
				make_factor <- TRUE
			}
		} else if (is.factor(values)) {
			levs <- levels(values)			
			values <- as.integer(values) - 1
			if (max(values, na.rm=TRUE) <= 255) {
				make_factor <- TRUE
			}
		}

		if (!(is.numeric(values) || is.integer(values) || is.logical(values))) {
			error("setValues", "values must be numeric, integer, logical, factor or character")
		}

		lv <- length(values)
		y <- rast(x)
		nc <- ncell(y)
		nl <- nlyr(y)
		opt <- spatOptions()
		
		if (lv == 1) {
			y@ptr$setValues(values, opt)
		} else {
			if (!((lv %% nc) == 0)) {
				warn("setValues", "length of values does not match the number of cells")
			}
			if (lv > (nc * nl)) {
				values <- values[1:(nc*nl)]
			} else if (lv < (nc * nl)) {
				values <- rep(values, length.out=nc*nl)
			}
			y@ptr$setValues(values, opt)
		}
		y <- messages(y)
		if (make_factor) {
			for (i in 1:nlyr(y)) {
				setCats(y, i, levs, 2)
			}
		}
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
	function(x) {
		src <- x@ptr$filenames
		src[src == ""] <= "memory"
		data.frame(source=src, nlyr=x@ptr$nlyrBySource(), stringsAsFactors=FALSE)
	}
)


setMethod("minmax", signature(x="SpatRaster"), 
	function(x) {
		have <- x@ptr$hasRange
		if (!any(have)) {
			warn("minmax", "min and max values not available. See 'setMinMax' or 'global'")
		}
		r <- rbind(x@ptr$range_min, x@ptr$range_max)
		r[,!have] <- c(Inf, -Inf)
		colnames(r) <- names(x)
		r
	}
)


setMethod("setMinMax", signature(x="SpatRaster"), 
	function(x, force=FALSE) {
		if (force) {
			x@ptr$setRange()
		} else if (any(!.hasMinMax(x))) {
			x@ptr$setRange()
		}
		x <- messages(x)
	}
)


setMethod("compareGeom", signature(x="SpatRaster", y="SpatRaster"), 
	function(x, y, ..., lyrs=FALSE, crs=TRUE, warncrs=FALSE, ext=TRUE, rowcol=TRUE, res=FALSE) {
		dots <- list(...)
		bool <- x@ptr$compare_geom(y@ptr, lyrs, crs, warncrs, ext, rowcol, res)
		messages(x, "compareGeom")
		if (length(dots)>1) {
			for (i in 1:length(dots)) {
				bool <- x@ptr$compare_geom(dots[[i]]@ptr, lyrs, crs, warncrs, ext, rowcol, res)
				messages(x, "compareGeom")
			}
		}
		bool
	}
)


setMethod("values", signature("SpatVector"), 
	function(x) {
		as.data.frame(x)
	}
)


setMethod("values<-", signature("SpatVector", "data.frame"), 
	function(x, value) {
		stopifnot(nrow(x) == nrow(value))
		x <- x[,0]
		# use cbind instead
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

setMethod("values<-", signature("SpatVector", "matrix"), 
	function(x, value) {
		`values<-`(x, data.frame(value))
	}
)


setMethod("values<-", signature("SpatVector", "ANY"),
	function(x, value) {
		if (!is.vector(value)) {
			error("values<-", "the values must be a data.frame, matrix or vector")
		}
		value <- rep(value, length.out=nrow(x))
		value <- data.frame(value=value)
		`values<-`(x, data.frame(value))
	}
) 



setMethod("values<-", signature("SpatVector", "NULL"), 
	function(x, value) {
		x@ptr$remove_df()
		x
	}
)

setMethod("setValues", signature("SpatVector", "ANY"), 
	function(x, values) {
		x@ptr <- x@ptr$deepcopy()
		`values<-`(x, values)
	}
)

