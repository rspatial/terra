# Author: Robert J. Hijmans
# Date :  June 2018
# Version 0.9
# License GPL v3

setMethod("hasValues", signature(x="SpatRaster"),
	function(x) {
		x@cpp$hasValues
	}
)


.makeDataFrame <- function(x, v, ...) {

	v <- data.frame(v, check.names=FALSE, ...)

#	factors=TRUE,
#	if (factors) {
	ff <- is.factor(x)
	if (any(ff)) {
		fs <- which(ff)
		cgs <- levels(x)
		for (f in fs) {
			cg <- cgs[[f]]
			i <- match(v[,f], cg[,1])
			if (!inherits(cg[[2]], "numeric")) {
				v[[f]] <- factor(cg[i, 2], levels=unique(cg[[2]]))
			} else {
				v[[f]] <- cg[i, 2]
			}
		}
	} 
	bb <- is.bool(x)
	if (any(bb)) {
		for (b in which(bb)) {
			v[[b]] = as.logical(v[[b]])
		}
	}
	ii <- (is.int(x) & (!ff) & (substr(datatype(x, TRUE), 1, 4) != "INT8"))
	if (any(ii)) {
		for (i in which(ii)) {
			v[[i]] = as.integer(v[[i]])
		}
	}
	dd <- !(bb | ii | ff)
	if (any(dd)) {
		d <- which(dd)
		v[,d] <- replace(v[,d], is.na(v[,d]), NA)
	}
	v
}


setMethod("readValues", signature(x="SpatRaster"),
function(x, row=1, nrows=nrow(x), col=1, ncols=ncol(x), mat=FALSE, dataframe=FALSE, ...) {
	stopifnot(row > 0 && nrows > 0)
	stopifnot(col > 0 && ncols > 0)
	v <- x@cpp$readValues(row-1, nrows, col-1, ncols)
	messages(x, "readValues")
	if (dataframe) {
		v <- matrix(v, ncol = nlyr(x))
		colnames(v) <- names(x)
		return(.makeDataFrame(x, v, ...) )
	} else if (mat) {
		if (all(is.int(x))) {
			v <- matrix(as.integer(v), ncol = nlyr(x))
		} else if (all(is.bool(x))) {
			v <- matrix(as.logical(v), ncol = nlyr(x))		
		} else {
			v <- matrix(v, ncol = nlyr(x))		
		}
		colnames(v) <- names(x)
	} else if (all(is.int(x))) {
		v <- as.integer(v)
	} else if (all(is.bool(x))) {
		v <- as.logical(v)
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
		setValues(x, value, keepnames=TRUE)
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
		m <- matrix(x@cpp$focalValues(w, fill, max(0, row-1), nrows, opt), ncol=prod(w), byrow=TRUE)
		messages(x, "focalValues")
		m
	}
)


mtrans <- function(mm, nc) {
	v <- NULL
	n <- ncol(mm) / nc
	for (i in 1:n) {
		j <- 1:nc + (i-1)*nc
		v <- c(v, as.vector(t(mm[, j])))
	}
	v
}


setMethod("setValues", signature("SpatRaster"),
	function(x, values, keeptime=TRUE, keepunits=TRUE, keepnames=FALSE, props=FALSE) {

		y <- rast(x, keeptime=keeptime, keepunits=keepunits, props=props)
		if (is.data.frame(values)) {
			# needs improvement to deal with mixed data types
			values <- as.matrix(values)
		}
		if (is.matrix(values)) {
			nl <- nlyr(x)
			d <- dim(values)
			if (!all(d == c(ncell(x), nl))) {
				ncx <- ncol(x)
				if ((d[1] == nrow(x)) && ((d[2] %% nl*ncx) == 0)) { 
					# raster-shaped matrix 
					if (ncx < d[2]) {
						values <- mtrans(values, ncx)
					} else {
						values <- as.vector(t(values))
					}
				} else if ((d[2] == nl) && (d[1] < ncell(x))) {
					if (d[1] > 1) warn("setValues", "values were recycled")
					values <- as.vector(apply(values, 2, function(i) rep_len(i, ncell(x))))
				} else {
					error("setValues","dimensions of the matrix do not match the SpatRaster")
				}
			} 
			if (!keepnames) {
				nms <- colnames(values)
				if (!is.null(nms)) names(y) <- nms
			}
		} else if (is.array(values)) {
			stopifnot(length(dim(values)) == 3)
			values <- as.vector(aperm(values, c(2,1,3)))
		}
		make_factor <- FALSE
		set_coltab <- FALSE
		if (is.character(values)) {
			if (all(substr(stats::na.omit(values), 1, 1) == "#")) {
				fv <- as.factor(values)
				if (length(levels(fv)) <= 256) {
					values <- as.integer(fv) #-1
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
				values <- as.integer(values) # -1 not needed anymore?
				make_factor <- TRUE
			}
		} else if (is.factor(values)) {
			levs <- levels(values)
			values <- as.integer(values)# - 1
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
			y@cpp$setValues(values, opt)
		} else {
			if (lv > (nc * nl)) {
				warn("setValues", "values is larger than the size of cells")
				values <- values[1:(nc*nl)]
			} else if (lv < (nc * nl)) {
				warn("setValues", "values were recycled")
				values <- rep(values, length.out=nc*nl)
			}
			y@cpp$setValues(values, opt)
		}
		y <- messages(y, "setValues")
		if (make_factor) {
			for (i in 1:nlyr(y)) {
				levs <- data.frame(value=1:length(levs), label=levs)
				set.cats(y, i, levs)
			}
			names(y) <- names(x)
		}
		if (set_coltab) {
			coltab(y) <- fv
		} else if (is.logical(values)) {
			if (!all(is.na(values))) y@cpp$setValueType(3)
		} else if (is.integer(values)) {
			y@cpp$setValueType(1)
		}
		y
	}
)



setMethod("inMemory", signature(x="SpatRaster"),
	function(x, bylayer=FALSE) {
		r <- x@cpp$inMemory
		if (bylayer) {
			nl <- .nlyrBySource(x)
			r <- rep(r, nl)
		}
		r
	}
)


#..hasValues <- function(x) { x@cpp$hasValues}
#..inMemory <- function(x) { x@cpp$inMemory }
#..filenames <- function(x) {	x@cpp$filenames }

subsetSource <- function(x, i) {
	x@cpp <- x@cpp$subsetSource(i-1)
	messages(x)
}

setMethod("sources", signature(x="SpatRaster"),
	function(x, nlyr=FALSE, bands=FALSE) {
		src <- x@cpp$filenames()
		Encoding(src) <- "UTF-8"
		if (bands) {
			nls <- x@cpp$nlyrBySource()
			d <- data.frame(sid=rep(1:length(src), nls),
						 source=rep(src, nls),
						 bands=x@cpp$getBands()+1, stringsAsFactors=FALSE)
			if (nlyr) {
				d$nlyr <- rep(nls, nls)
			}
			d
		} else if (nlyr) {
			data.frame(source=src, nlyr=x@cpp$nlyrBySource(), stringsAsFactors=FALSE)
		} else {
			src
		}
	}
)

setMethod("sources", signature(x="SpatRasterCollection"),
	function(x, nlyr=FALSE, bands=FALSE) {
		if (nlyr | bands) {
			x <- lapply(x, function(i) sources(i, nlyr, bands))
			x <- lapply(1:length(x), function(i) cbind(cid=i, x[[i]]))
			do.call(rbind, x)
		} else {
			sapply(x, sources)
		}
	}
)

setMethod("sources", signature(x="SpatRasterDataset"),
	function(x, nlyr=FALSE, bands=FALSE) {
		if (nlyr | bands) {
			x <- lapply(x, function(i) sources(i, nlyr, bands))
			x <- lapply(1:length(x), function(i) cbind(cid=i, x[[i]]))
			do.call(rbind, x)
		} else {
			x@cpp$filenames()
		}
	}
)


setMethod("sources", signature(x="SpatVector"),
	function(x) {
		if (x@cpp$source != "") {
			if (x@cpp$layer != tools::file_path_sans_ext(basename(x@cpp$source))) {
				paste0(x@cpp$source, "::", x@cpp$layer)
			} else {
				x@cpp$source
			}
		} else {
			""
		}
	}
)

setMethod("sources", signature(x="SpatVectorProxy"),
	function(x) {
		if (x@cpp$v$layer != tools::file_path_sans_ext(basename(x@cpp$v$source))) {
			paste0(x@cpp$v$source, "::", x@cpp$v$layer)
		} else {
			x@cpp$v$source
		}
	}
)

setMethod("hasMinMax", signature(x="SpatRaster"),
	function(x) {
		x@cpp$hasRange
	}
)

setMethod("minmax", signature(x="SpatRaster"),
	function(x, compute=FALSE) {
		have <- x@cpp$hasRange
		if (!all(have)) {
			if (compute) {
				opt <- spatOptions()
				x@cpp$setRange(opt, FALSE)
			} else {
				warn("minmax", "min and max values not available for all layers. See 'setMinMax' or 'global'")
			}
		}
		r <- rbind(x@cpp$range_min, x@cpp$range_max)
		if (!compute) {
			r[,!have] <- c(Inf, -Inf)
		}
		colnames(r) <- names(x)
		rownames(r) <- c("min", "max")
		r
	}
)


setMethod("setMinMax", signature(x="SpatRaster"),
	function(x, force=FALSE) {
		opt <- spatOptions()
		if (force) {
			x@cpp$setRange(opt, TRUE)
		} else if (!all(hasMinMax(x))) {
			x@cpp$setRange(opt, FALSE)
		}
		x <- messages(x, "setMinMax")
	}
)



setMethod("compareGeom", signature(x="SpatRaster", y="SpatRaster"),
	function(x, y, ..., lyrs=FALSE, crs=TRUE, warncrs=FALSE, ext=TRUE, rowcol=TRUE, res=FALSE, stopOnError=TRUE, messages=FALSE) {
		opt <- spatOptions("")
		out <- x@cpp$compare_geom(y@cpp, lyrs, crs, opt$tolerance, warncrs, ext, rowcol, res)
		if (stopOnError) {
			messages(x, "compareGeom")
		} else {
			m <- NULL
			if (x@cpp$has_warning()) {
				m <- x@cpp$getWarnings()
			}
			if (x@cpp$has_error()) {
				m <- c(m, x@cpp$getError())
			}
			if (!is.null(m) && messages) {
				message(paste(m, collapse="\n"))
			}
		}
		dots <- list(...)
		if (length(dots) > 0) {
			for (i in 1:length(dots)) {
				if (!inherits(dots[[i]], "SpatRaster")) {
					error("compareGeom", "all additional arguments must be a SpatRaster")
				}
				out <- out & 
					compareGeom(x, dots[[i]], lyrs=lyrs, crs=crs, warncrs=warncrs, ext=ext, rowcol=rowcol, res=res, stopOnError=stopOnError, messages=messages)
			}
		}
		out
	}
)


setMethod("compareGeom", signature(x="SpatRaster", y="list"),
	function(x, y, ..., lyrs=FALSE, crs=TRUE, warncrs=FALSE, ext=TRUE, rowcol=TRUE, res=FALSE, stopOnError=TRUE, messages=FALSE) {
		dots <- list(...)
		if (length(dots) > 0) {
			y <- c(y, dots)
		}
		isr <- sapply(y, inherits, "SpatRaster")
		if (!all(isr)) {
			y <- y[isr]
			if (length(y) == 0) error("compareGeom", "none of the elements of y is a SpatRaster")
			n <- sum(!isr)
			if (n > 1) {
				warn("compareGeom", paste(n, "elements of y are not a SpatRaster"))
			} else {
				warn("compareGeom", paste("1 element of y is not a SpatRaster"))
			}
		}
		out <- sapply(y, compareGeom, y=x, lyrs=lyrs, crs=crs, warncrs=warncrs, ext=ext, rowcol=rowcol, res=res, stopOnError=stopOnError, messages=messages)
		if (!all(isr)) {
			res <- rep(NA, length(isr))
			res[isr] <- out
			return(res)
		}
		all(out)
	}
)


setMethod("compareGeom", signature(x="SpatRasterCollection", y="missing"),
	function(x, y, ..., lyrs=FALSE, crs=TRUE, warncrs=FALSE, ext=TRUE, rowcol=TRUE, res=FALSE, stopOnError=FALSE, messages=FALSE) {
		x <- as.list(x)
		dots <- list(...)
		if (length(dots) > 0) {
			error("compareGeom", "when x is a SpatRasterCollection, additional arguments are ignored")
		}
		if (length(x) == 1) return(TRUE)
		out <- sapply(x[-1], compareGeom, y=x[1], lyrs=lyrs, crs=crs, warncrs=warncrs, ext=ext, rowcol=rowcol, res=res, stopOnError=stopOnError, messages=messages)
		out
	}
)


setMethod("compareGeom", signature(x="SpatRaster", y="SpatRasterCollection"),
	function(x, y, ..., lyrs=FALSE, crs=TRUE, warncrs=FALSE, ext=TRUE, rowcol=TRUE, res=FALSE, stopOnError=FALSE, messages=FALSE) {
		dots <- list(...)
		if (length(dots) > 0) {
			error("compareGeom", "when y is a SpatRasterCollection, additional arguments are ignored")
		}
		out <- sapply(as.list(y), compareGeom, y=x, lyrs=lyrs, crs=crs, warncrs=warncrs, ext=ext, rowcol=rowcol, res=res, stopOnError=stopOnError, messages=messages)
		all(out)
	}
)


setMethod("compareGeom", signature(x="SpatVector", y="SpatVector"),
	function(x, y, tolerance=0) {
		out <- x@cpp$equals_between(y@cpp, tolerance)
		x <- messages(x, "compareGeom")
		out[out == 2] <- NA
		matrix(as.logical(out), nrow=nrow(x), byrow=TRUE)
	}
)

setMethod("compareGeom", signature(x="SpatVector", y="SpatVector"),
	function(x, y, tolerance=0) {
		out <- x@cpp$equals_within(tolerance)
		x <- messages(x, "compareGeom")
		out[out == 2] <- NA
		out <- matrix(as.logical(out), nrow=nrow(x), byrow=TRUE)
		out
	}
)


setMethod("all.equal", signature(target="SpatRaster", current="SpatRaster"),
	function(target, current, maxcell=10000, ...) {
		a <- base::all.equal(rast(target), rast(current))
		if (isTRUE(a) && maxcell > 0) {
			hvT <- hasValues(target)
			hvC <- hasValues(current)
			if (hvT && hvC) {
				s1 <- spatSample(current, maxcell, "regular")
				s2 <- spatSample(target, maxcell, "regular")
				a <- all.equal(s1, s2, ...)
			} else if (hvT || hvC) {
				if (hvT) {
					a <- "target has cell values, current does not"
				} else {
					a <- "current has cell values, target does not"				
				}
			}
		}
		a
	}
)



setMethod("values", signature("SpatVector"),
	function(x, ...) {
		as.data.frame(x, ...)
	}
)


setMethod("values<-", signature("SpatVector", "data.frame"),
	function(x, value) {
		x@cpp <- x@cpp$deepcopy()
		if (ncol(value) == 0) {
			x@cpp$remove_df()
			return(x)
		}
		value <- .makeSpatDF(value)
		x@cpp$set_df(value)
		messages(x)
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
		x@cpp <- x@cpp$deepcopy()
		x@cpp$remove_df()
		x
	}
)

setMethod("setValues", signature("SpatVector"),
	function(x, values) {
		x@cpp <- x@cpp$deepcopy()
		`values<-`(x, values)
	}
)

