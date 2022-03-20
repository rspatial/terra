# Author: Robert J. Hijmans
# Date :  June 2008
# Version 0.9
# License GPL v3

setMethod("hasValues", signature(x="SpatRaster"), 
	function(x) {
		x@ptr$hasValues
	}
)


.makeDataFrame <- function(x, v, factors=TRUE, ...) {

	v <- data.frame(v, check.names=FALSE, ...)

	if (factors) {
		ff <- is.factor(x)
		if (any(ff)) {
			ff <- which(ff)
			cgs <- cats(x)
			for (f in ff) {
				cg <- cgs[[f]]
				i <- match(v[,f], cg[,1])
				act <- activeCat(x, f) + 1
				if (!inherits(cg[[act]], "numeric")) {
					v[[f]] <- factor(cg[i, act], levels=unique(cg[[act]]))
				} else {
					v[[f]] <- cg[i, act]
				}
			}
		} else {
			bb <- is.bool(x)
			if (any(bb)) {
				for (b in which(bb)) {
					v[[b]] = as.logical(v[[b]])
				}
			}

			ii <- is.int(x)
			if (any(ii)) {
				for (i in which(ii)) {
					v[[i]] = as.integer(v[[i]])
				}
			}
		}		
	} else {
		bb <- is.bool(x)
		if (any(bb)) {
			for (b in which(bb)) {
				v[[b]] = as.logical(v[[b]])
			}
		}

		ii <- is.int(x)
		if (any(ii)) {
			for (i in which(ii)) {
				v[[i]] = as.integer(v[[i]])
			}
		}
	}
	v
}


setMethod("readValues", signature(x="SpatRaster"), 
function(x, row=1, nrows=nrow(x), col=1, ncols=ncol(x), mat=FALSE, dataframe=FALSE, ...) {
	stopifnot(row > 0 && nrows > 0)
	stopifnot(col > 0 && ncols > 0)
	v <- x@ptr$readValues(row-1, nrows, col-1, ncols)
	messages(x, "readValues")
	if (dataframe || mat) {
		v <- matrix(v, ncol = nlyr(x))
		colnames(v) <- names(x)
		if (dataframe) {
			return(.makeDataFrame(x, v, factors=TRUE, ...) )
		}
	}
	v
}
)


setMethod("values", signature(x="SpatRaster"), 
function(x, mat=TRUE, dataframe=FALSE, row=1, nrows=nrow(x), col=1, ncols=ncol(x), na.rm=FALSE, ...) {
	readStart(x)
	on.exit(readStop(x))
	v <- readValues(x, row, nrows, col, ncols, mat=mat, dataframe=dataframe, ...)
	messages(x, "values")
	if (na.rm) {
		if (!is.null(dim(v))) {
			v[stats::complete.cases(v), , drop=FALSE]
		} else {
			v[!is.na(v)]
		}
	} else {
		v
	}
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
		opt <- spatOptions()
		m <- matrix(x@ptr$focalValues(w, fill, max(0, row-1), nrows, opt), ncol=prod(w), byrow=TRUE)
		messages(x, "focalValues")
		m
	}
)


setMethod("setValues", signature("SpatRaster"), 
	function(x, values, keeptime=TRUE, keepunits=TRUE, props=FALSE) {

		y <- rast(x, keeptime=keeptime, keepunits=keepunits, props=props)

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
		set_coltab <- FALSE
		if (is.character(values)) {
			if (all(substr(na.omit(values), 1, 1) == "#")) {
				fv <- as.factor(values)
				if (length(levels(fv)) <= 256) {
					values <- as.integer(fv)-1
					fv <- levels(fv)
					set_coltab <- TRUE
				} else {
					fv <- NULL
					values <- t(grDevices::col2rgb(values))
					y <- rast(y, nlyr=3, names=c("red", "green", "blue"))
					RGB(y) <- 1:3
				}
			} else {
				values <- as.factor(values)
				levs <- levels(values)
				values <- as.integer(values) - 1 # -1 not needed anymore?
				make_factor <- TRUE
			}
		} else if (is.factor(values)) {
			levs <- levels(values)
			values <- as.integer(values) - 1
			make_factor <- TRUE
		} 

		if (!(is.numeric(values) || is.integer(values) || is.logical(values))) {
			error("setValues", "values must be numeric, integer, logical, factor or character")
		}

		lv <- length(values)
		nc <- ncell(y)
		nl <- nlyr(y)
		opt <- spatOptions()

		if (lv == 1) {
			y@ptr$setValues(values, opt)
		} else {
			if (lv > (nc * nl)) {
				warn("setValues", "values is larger than the size of cells")
				values <- values[1:(nc*nl)]
			} else if (lv < (nc * nl)) {
				warn("setValues", "values were recycled")
				values <- rep(values, length.out=nc*nl)
			} 
			y@ptr$setValues(values, opt)
		}
		y <- messages(y, "setValues")
		if (make_factor) {
			for (i in 1:nlyr(y)) {
				set.cats(y, i, levs, 2)
			}
		}
		if (set_coltab) {
			coltab(y) <- fv
		} else if (is.logical(values)) {
			y@ptr$setValueType(3)
		} else if (is.integer(values)) {
			y@ptr$setValueType(1)
		}
		y
	}
)



setMethod("inMemory", signature(x="SpatRaster"), 
	function(x, bylayer=FALSE) {
		r <- x@ptr$inMemory
		if (bylayer) {
			nl <- .nlyrBySource(x)
			r <- rep(r, nl)
		}
		r
	}
)


#..hasValues <- function(x) { x@ptr$hasValues}
#..inMemory <- function(x) { x@ptr$inMemory }
#..filenames <- function(x) {	x@ptr$filenames }


setMethod("sources", signature(x="SpatRaster"), 
	function(x, nlyr=FALSE, bands=FALSE) {
		src <- x@ptr$filenames
		Encoding(src) <- "UTF-8"
		if (bands) {
			nls <- x@ptr$nlyrBySource()
			d <- data.frame(sid=rep(1:length(src), nls), 
						 source=rep(src, nls), 
						 bands=x@ptr$getBands()+1, stringsAsFactors=FALSE)
			if (nlyr) {
				d$nlyr <- rep(nls, nls) 
			}
			d			 
		} else if (nlyr) {
			data.frame(source=src, nlyr=x@ptr$nlyrBySource(), stringsAsFactors=FALSE)
		} else {
			src
		}
	}
)

setMethod("hasMinMax", signature(x="SpatRaster"), 
	function(x) {
		x@ptr$hasRange
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
		opt <- spatOptions()
		if (force) {
			x@ptr$setRange(opt)
		} else if (!all(hasMinMax(x))) {
			x@ptr$setRange(opt)
		}
		x <- messages(x, "setMinMax")
	}
)


setMethod("compareGeom", signature(x="SpatRaster", y="SpatRaster"), 
	function(x, y, ..., lyrs=FALSE, crs=TRUE, warncrs=FALSE, ext=TRUE, rowcol=TRUE, res=FALSE, stopOnError=TRUE, messages=TRUE) {
		dots <- list(...)
		opt <- spatOptions("")
		res <- x@ptr$compare_geom(y@ptr, lyrs, crs, opt$tolerance, warncrs, ext, rowcol, res)
		if (stopOnError) {
			messages(x, "compareGeom")
		} else {
			m <- NULL
			if (x@ptr$has_warning()) { 
				m <- x@ptr$getWarnings()
			}
			if (x@ptr$has_error()) {
				m <- c(m, x@ptr$getError())
			}
			if (!is.null(m) && messages) {
				message(paste(m, collapse="\n"))
			}
		}
		if (length(dots)>1) {
			for (i in 1:length(dots)) {
				bool <- x@ptr$compare_geom(dots[[i]]@ptr, lyrs, crs, warncrs, ext, rowcol, res)
				if (stopOnError) messages(x, "compareGeom")
				res <- bool & res
			}
		}
		res
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

setMethod("setValues", signature("SpatVector"), 
	function(x, values) {
		x@ptr <- x@ptr$deepcopy()
		`values<-`(x, values)
	}
)

