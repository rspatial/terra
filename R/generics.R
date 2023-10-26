# Author: Robert J. Hijmans
# Date : October 2018
# Version 1.0
# License GPL v3


setMethod("all.equal", signature(target="SpatRaster", current="SpatRaster"),
	function(target, current, maxcell=10000, ...) {
		first <- all.equal.default(target, current)
		if (isTRUE(first)) {
			if (hasValues(target)) {
				if (ncell(x) > maxcell) {
					s <- round(100 * maxcell / ncell(x))
					warn("all.equal", paste0("using a sample of ", s, "% of the cells"))
				}
				vt <- spatSample(target, maxcell, "regular")
				ct <- spatSample(current, maxcell, "regular")
				all.equal(vt, ct, ...)
			} else {
				first
			}
		} else {
			first
		}
	}
)



setMethod("weighted.mean", signature(x="SpatRaster", w="numeric"),
	function(x, w, na.rm=FALSE, filename="", ...) {
		opt <- spatOptions(filename, ...)
		x@ptr <- x@ptr$wmean_vect(w, na.rm, opt)
		messages(x, "weighted.mean")
	}
)


setMethod("weighted.mean", signature(x="SpatRaster", w="SpatRaster"),
	function(x, w, na.rm=FALSE, filename="", ...) {
		opt <- spatOptions(filename, ...)
		x@ptr <-x@ptr$wmean_rast(w@ptr, na.rm, opt)
		messages(x, "weighted.mean")
	}
)



setMethod("patches", signature(x="SpatRaster"),
	function(x, directions=4, zeroAsNA=FALSE, allowGaps=TRUE, filename="", ...) {
		if (allowGaps) {
			opt <- spatOptions(filename, ...)
			x@ptr <- x@ptr$patches(directions[1], zeroAsNA[1], opt)
			messages(x, "patches")
		} else {
			opt <- spatOptions()
			x@ptr <- x@ptr$patches(directions[1], zeroAsNA[1], opt)
			x <- messages(x, "patches")
			f <- freq(x)
			fr <- cbind(f[,2], 1:nrow(f))
			fr <- fr[fr[,1] != fr[,2], ,drop=FALSE]
			if (nrow(fr) > 0) {
				x <- classify(x, fr, filename=filename, ...)
			} else if (filename != "") {
				x <- writeRaster(x, filename=filename, ...)
			}
			x
		}
	}
)


setMethod("origin", signature(x="SpatRaster"),
	function(x) {
		x@ptr$origin
	}
)


setMethod("origin<-", signature("SpatRaster"),
	function(x, value) {
		value <- rep(value, length.out=2)
		dif <- value - origin(x)
		res <- res(x)
		dif[1] <- dif[1] %% res[1]
		dif[2] <- dif[2] %% res[2]
		for (i in 1:2) {
			if (dif[i] < 0) {
				if ((dif[i] + res[i]) < abs(dif[i])) {
					dif[i] <- dif[i] + res[i]
				}
			} else {
				if (abs(dif[i] - res[i]) < dif[i]) {
					dif[i] <- dif[i] - res[i]
				}
			}
		}
		e <- as.vector(ext(x))
		e["xmin"] <- e["xmin"] + dif[1]
		e["xmax"] <- e["xmax"] + dif[1]
		e["ymin"] <- e["ymin"] + dif[2]
		e["ymax"] <- e["ymax"] + dif[2]
		ext(x) <- e
		return(x)
	}
)




setMethod("align", signature(x="SpatExtent", y="SpatRaster"),
	function(x, y, snap="near") {
		x@ptr <- y@ptr$align(x@ptr, tolower(snap))
		#messages(x, "align")
		x
	}
)

setMethod("align", signature(x="SpatExtent", y="numeric"),
	function(x, y) {
		x@ptr <- x@ptr$align(y, "")
		x
	}
)

setMethod("cellSize", signature(x="SpatRaster"),
	function(x, mask=TRUE, unit="m", transform=TRUE, rcx=100, filename="", ...) {
		opt <- spatOptions(filename, ...)
		x@ptr <- x@ptr$rst_area(mask, unit, transform, rcx, opt)
		messages(x, "cellSize")
	}
)



setMethod ("expanse", "SpatRaster",
	function(x, unit="m", transform=TRUE) {

		byvalue = FALSE
		opt <- spatOptions()
		if (byvalue) {
			v <- x@ptr$area_by_value(opt)
			x <- messages(x, "expanse")
			v <- lapply(1:length(v), function(i) cbind(i, matrix(v[[i]], ncol=2)))
			v <- do.call(rbind, v)
			colnames(v) <- c("layer", "value", "area")
		} else {
			v <- x@ptr$sum_area(unit, transform, opt)
			x <- messages(x, "expanse")
		}
		return(v)
	}
)



setMethod("atan2", signature(y="SpatRaster", x="SpatRaster"),
	function(y, x) {
		opt <- spatOptions(filename="", overwrite=TRUE)
		y@ptr <- y@ptr$atan2(x@ptr, opt)
		messages(y, "atan2")
	}
)

setMethod("atan_2", signature(y="SpatRaster", x="SpatRaster"),
	function(y, x, filename="", ...) {
		opt <- spatOptions(filename=filename, ...)
		y@ptr <- y@ptr$atan2(x@ptr, opt)
		messages(y, "atan_2")
	}
)


setMethod("boundaries", signature(x="SpatRaster"),
	function(x, classes=FALSE, inner=TRUE, directions=8, falseval=0, filename="", ...) {
		opt <- spatOptions(filename, ...)
		type <- ifelse(inner[1], "inner", "outer")
		x@ptr <- x@ptr$boundaries(classes[1], type, directions[1], falseval[1], opt)
		messages(x, "boundaries")
	}
)


.collapseSources <- function(x) {
	x@ptr <- x@ptr$collapse_sources()
	messages(x, "tighten")
}

setMethod("deepcopy", signature("SpatRaster"),
	function(x) {
		x@ptr <- x@ptr$deepcopy()
		x
	}
)



setMethod("split", signature(x="SpatRaster"),
	function(x, f) {
		stopifnot(length(f) == nlyr(x))
		stopifnot(!any(is.na(f)))
		u <- unique(f)
		lapply(u, function(i) x[[f==i]])
	}
)


setMethod("add<-", signature("SpatRaster", "SpatRaster"),
	function(x, value) {
		if (x@ptr$same(value@ptr)) {
			x@ptr <- x@ptr$deepcopy()
		}
		opt <- spatOptions()
		x@ptr$addSource(value@ptr, FALSE, opt)
		messages(x, "add")
	}
)

setMethod("tighten", signature("SpatRaster"),
	function(x) {
		x@ptr <- x@ptr$collapse_sources()
		messages(x, "tighten")
	}
)

setMethod("tighten", signature("SpatRasterDataset"),
	function(x) {
		y <- new("SpatRaster")
		y@ptr <- x@ptr$collapse()
		messages(y, "tighten")
	}
)


#setMethod("c", signature(x="SpatRaster"),
#	function(x, ...) {
#		s <- sds(list(x, ...))
#		x@ptr <- s@ptr$collapse()
#		x <- messages(x, "c")
#		try( x@ptr <- x@ptr$collapse_sources() )
#		messages(x, "c")
#	}
#)



#cbind.SpatVector <- function(x, y, ...) {
#	if (inherits(y, "SpatVector")) {
#		y <- y@ptr$df
#	} else {
#		stopifnot(inherits(y, "data.frame"))
#		y <- terra:::.makeSpatDF(y)
#	}
#	x@ptr <- x@ptr$cbind(y)
#	messages(x, "cbind")
#}

cbind.SpatVector <- function(x, y, ...) {
	dots <- list(y, ...)
	for (y in dots) {
		if (inherits(y, "SpatVector")) {
			y <- y@ptr$df
		} else {
			stopifnot(inherits(y, "data.frame"))
			y <- .makeSpatDF(y)
		}
		x@ptr <- x@ptr$cbind(y)
		x <- messages(x, "cbind")
	}
	x
}

rbind.SpatVector <- function(x, y, ...) {
	skipped <- FALSE
	if (missing(y) || is.null(y)) return(x)
	stopifnot(inherits(y, "SpatVector"))
	s <- svc(x, y, ...)
	x@ptr <- s@ptr$append()
	messages(x, "rbind")

#	dots <- list(...)
#	if (is.null(dots)) {
#		x@ptr <- x@ptr$rbind(y@ptr, FALSE)
#	} else {

	#if (!is.null(dots)) {
	#	for (y in dots) {
	#		stopifnot(inherits(y, "SpatVector"))
	#		x@ptr <- x@ptr$rbind(y@ptr, FALSE)
	#		x <- messages(x, "rbind")
	#	}
	#}
	#x
}


# this way names of named arguments are used
setMethod("c", signature(x="SpatRaster"),
	function(x, ..., warn=TRUE) {
		if (missing(x)) {
			rast(list(...), warn=warn)
		} else {
			rast(list(x, ...), warn=warn)
		}
	}
)

#setMethod("c", signature(x="SpatRaster"),
#	function(x, ..., warn=TRUE) {
#		skips <- 0
#		hv <- hasValues(x)
#		dots <- list(...)
#		x@ptr <- x@ptr$deepcopy()
#		opt <- spatOptions()
#		for (i in dots) {
#			if (inherits(i, "SpatRaster")) {
#				x@ptr$addSource(i@ptr, warn, opt)
#				if (x@ptr$messages$has_error) {
#					messages(x, "c")
#					return()
#				}
#			} else {
#				skips = skips + 1
#			}
#		}
#		if (skips > 0) warn("c,SpatRaster", paste("skipped", skips, "object(s) that are not SpatRaster"))
#		messages(x, "c")
#		x
#	}
#)



setMethod("rep", signature(x="SpatRaster"),
	function(x, ...) {
		i <- rep(1:nlyr(x), ...)
		x[[i]]
	}
)


setMethod("clamp", signature(x="SpatRaster"),
	function(x, lower=-Inf, upper=Inf, values=TRUE, filename="", ...) {
		opt <- spatOptions(filename, ...)
		x@ptr <- x@ptr$clamp(lower, upper, values[1], opt)
		messages(x, "clamp")
	}
)


setMethod("clamp", signature(x="numeric"),
function(x, lower=-Inf, upper=Inf, values=TRUE, ...) {
	stopifnot(lower <= upper)
	if (values) {
		x[x < lower] <- lower
		x[x > upper] <- upper
	} else {
		x[x < lower] <- NA
		x[x > upper] <- NA
	}
	x
}
)

setMethod("classify", signature(x="SpatRaster"),
function(x, rcl, include.lowest=FALSE, right=TRUE, others=NULL, brackets=TRUE, filename="", ...) {

	bylayer = FALSE

	if (is.data.frame(rcl)) {
		rcl <- as.matrix(rcl)
	}

	right <- ifelse(is.na(right), 2, ifelse(isTRUE(right), 1, 0))
	include.lowest <- as.logical(include.lowest[1])

	opt <- spatOptions(filename, ...)
	if (is.null(others)) {
		others <- FALSE
		othersValue <- 0
	} else {
		othersValue <- others[1]
		others <- TRUE
	}
	keepcats <- FALSE
	if (inherits(rcl[1], "character")) {
		if (nlyr(x) > 1) {
			error("classify", "rcl has characters. That is not allowed with multiple layers in x")
		}
		if (!is.factor(x)) {
			error("classify", "rcl has characters but x is not categorical")
		}
		if (ncol(rcl) != 2) {
			error("classify", "rcl has characters. It should have 2 columns")
		}
		levs <- cats(x, 1, active=TRUE)[[1]]
		rc1 <- levs[,1][match(rcl[,1], levs[,2])]
		rc2 <- levs[,1][match(rcl[,2], levs[,2])]
		rcl <- cbind(rc1, rc2)[!is.na(rc1), ]
		keepcats <- TRUE
	}

    x@ptr <- x@ptr$classify(as.vector(rcl), NCOL(rcl), right, include.lowest, others, othersValue, bylayer[1], brackets[1], keepcats, opt)
	messages(x, "classify")
}
)

setMethod("subst", signature(x="SpatRaster"),
function(x, from, to, filename="", ...) {
	opt <- spatOptions(filename, ...)

	if (inherits(from, "data.frame")) {
		from <- as.matrix(from)
	}
	if (inherits(to, "data.frame")) {
		to <- as.matrix(to)
	}
	tom <- inherits(to, "matrix")
	frm <- inherits(from, "matrix")

	keepcats <- FALSE
	fromc <- inherits(from[1], "character")
	toc <- inherits(to[1], "character")
	if (fromc || toc) {
		if (!(fromc && toc)) {
			error("subst", "either both or neither from and to should have character values")
		}
		if (!is.factor(x)) {
			error("subst", "from has characters but x is not categorical")
		}
		if (nlyr(x) > 1) {
			error("subst", "you can only use characters if x has 1 layer")
		}
		if (inherits(to, "matrix")) {
			if (ncol(to) == 1) {
				to <- as.vector(to)
			} else if (ncol(to) != nlyr(x)) {
				to <- as.vector(to[,1])
				warn("subst", "only the first column of 'to' is used with factors")
			}
		}
		levs <- cats(x, 1, active=TRUE)[[1]]
		from <- levs[,1][match(from, levs[,2])]
		to <- levs[,1][match(to, levs[,2])]
		keepcats <- TRUE
	}

	if (tom && frm) {
		error("subst", "either 'to' or 'from' can be a matrix, not both")
	}

	if (tom) {
		nms <- colnames(to)
		if (!is.null(nms)) 	opt$names = nms
		x@ptr <- x@ptr$replaceValues(from, to, ncol(to), keepcats, opt)
	} else if (frm) {
		x@ptr <- x@ptr$replaceValues(as.vector(t(from)), to, -ncol(from), FALSE, opt)
	} else {
		x@ptr <- x@ptr$replaceValues(from, to, 0, keepcats, opt)
	}
	messages(x, "subst")
}
)


.getExt <- function(y, method="crop") {
	if (!inherits(y, "SpatExtent")) {
		e <- try(ext(y), silent=TRUE)
		if (inherits(e, "try-error")) {
			e <- try(ext(vect(y)), silent=TRUE)
			if (inherits(e, "try-error")) {
				error(method, "cannot get a SpatExtent from y")
			}
		}
		y <- e
	}
	y
}

setMethod("crop", signature(x="SpatRaster", y="ANY"),
	function(x, y, snap="near", mask=FALSE, touches=TRUE, filename="", ...) {
		opt <- spatOptions(filename, ...)
		if (mask && inherits(y, "SpatVector")) {
			x@ptr <- x@ptr$crop_mask(y@ptr, snap[1], touches[1], opt)
		} else {
			y <- .getExt(y, method="crop")
			x@ptr <- x@ptr$crop(y@ptr, snap[1], opt)
		}
		messages(x, "crop")
	}
)

setMethod("crop", signature(x="SpatRasterDataset", y="ANY"),
	function(x, y, snap="near", filename="", ...) {
		if (all(filename != "")) {
			ext <- tools::file_ext(filename)
			filename <- tools::file_path_sans_ext(filename)
			filename <- paste0(make.unique(filename, sep="_"), ext)
		}
		opt <- spatOptions(filename, ...)
		y <- .getExt(y, method="crop")
		x@ptr <- x@ptr$crop(y@ptr, snap[1], opt)
		messages(x, "crop")
	}
)



setMethod("selectRange", signature(x="SpatRaster"),
	function(x, y, z=1, repint=0, filename="", ...) {
		opt <- spatOptions(filename, ...)
		x@ptr <- x@ptr$selRange(y@ptr, z, repint, opt)
		messages(x, "selectRange")
	}
)

setMethod("cover", signature(x="SpatRaster", y="SpatRaster"),
	function(x, y, values=NA, filename="", ...) {
		opt <- spatOptions(filename, ...)
		x@ptr <- x@ptr$cover(y@ptr, values, opt)
		messages(x, "cover")
	}
)


setMethod("diff", signature(x="SpatRaster"),
	function(x, lag=1, filename="", ...) {
		n <- nlyr(x)
		lag <- round(lag)
		if ((lag < 1) | (lag >= n)) {
			error("diff", "lag must be > 0 and < nlyr(x)")
		}
		y <- x[[-((n-lag+1):n)]]
		x <- x[[-(1:lag)]]
		opt <- spatOptions(filename, ...)
		x@ptr <- x@ptr$arith_rast(y@ptr, "-", opt)
		messages(x, "diff")
	}
)


setMethod("disagg", signature(x="SpatRaster"),
	function(x, fact, method="near", filename="", ...) {
		stopifnot(method %in% c("near", "bilinear"))
		if (method == "bilinear") {
			y <- disagg(rast(x), fact)
			r <- resample(x, y, "bilinear", filename=filename, ...)
			return(r)
		}
		opt <- spatOptions(filename, ...)
		x@ptr <- x@ptr$disaggregate(fact, opt)
		messages(x, "disagg")
	}
)


setMethod("flip", signature(x="SpatRaster"),
	function(x, direction="vertical", filename="", ...) {
		d <- match.arg(direction, c("vertical", "horizontal"))
		opt <- spatOptions(filename, ...)
		x@ptr <- x@ptr$flip(d == "vertical", opt)
		messages(x, "flip")
	}
)


setMethod("freq", signature(x="SpatRaster"),
	function(x, digits=0, value=NULL, bylayer=TRUE, usenames=FALSE) {

		opt <- spatOptions()
		if (!bylayer) usenames <- FALSE

		if (!is.null(value)) {
			value <- unique(value)
			if (length(value) > 1) {
				error("freq", "value must have a length of one")
			}
			if (is.character(value)) {
				value <- value[value != ""]
				if (length(value) == 0) {
					error("freq", "no valid value")
				}
				ff <- is.factor(x)
				if (!any(ff)) {
					error("freq", "a character value is only meaningful for categorical rasters")
				}
				f <- freq(x[[ff]])
				if (usenames) {
					f$layer <- names(x)[f$layer]
				}
				f <- f[f$label == value,]
				return(f)
			}

			if (is.na(digits)) {
				v <- x@ptr$count(value, bylayer[1], FALSE, 0, opt)
			} else {
				v <- x@ptr$count(value, bylayer[1], TRUE, digits, opt)
				value <- round(value, digits)
			}
			if (bylayer) {
				v <- data.frame(layer=1:nlyr(x), value=value, count=v)
			} else {
				v <- data.frame(value=value, count=v)
			}

		} else {
			if (is.na(digits)) {
				v <- x@ptr$freq(bylayer[1], FALSE, 0, opt)
			} else {
				v <- x@ptr$freq(bylayer[1], TRUE, digits, opt)
			}
			v <- lapply(1:length(v), function(i) cbind(i, matrix(v[[i]], ncol=2)))
			v <- do.call(rbind, v)
			v <- as.data.frame(v)
			colnames(v) <- c("layer", "value", "count")
			ff <- is.factor(x)
			if (any(ff)) {
				cgs <- cats(x)
				v <- data.frame(v)
				for (f in which(ff)) {
					cg <- cgs[[f]]
					j <- which(v[,1] == f)
					i <- match(v[j,2], cg[,1])
					act <- activeCat(x, f) + 1
					if (!inherits(cg[[act]], "numeric")) {
						v[j, 2] <- as.character(factor(cg[i, act], levels=unique(cg[[act]])))
					} else {
						v[j, 2] <- cg[i, act]
					}
				}
			}
			if (nlyr(x) > 1 && !bylayer) {
				v <- aggregate(v[,"count",drop=FALSE], v[,"value", drop=FALSE], sum)
			}
		}
		if (usenames) {
			v$layer <- names(x)[v$layer]
		}
		v
	}
)



setMethod("mask", signature(x="SpatRaster", mask="SpatRaster"),
	function(x, mask, inverse=FALSE, maskvalues=NA, updatevalue=NA, filename="", ...) {
		opt <- spatOptions(filename, ...)
		x@ptr <- x@ptr$mask_raster(mask@ptr, inverse[1], maskvalues, updatevalue[1], opt)
		messages(x, "mask")
	}
)

setMethod("mask", signature(x="SpatRaster", mask="SpatVector"),
	function(x, mask, inverse=FALSE, updatevalue=NA, touches=TRUE, filename="", ...) {
		opt <- spatOptions(filename, ...)
		x@ptr <- x@ptr$mask_vector(mask@ptr, inverse[1], updatevalue[1], touches[1], opt)
		messages(x, "mask")
	}
)

setMethod("mask", signature(x="SpatRaster", mask="sf"),
   function(x, mask, inverse=FALSE, updatevalue=NA, touches=TRUE, filename="", ...) {
		mask <- vect(mask)
		mask(x, mask, inverse=inverse, updatevalue=updatevalue, touches=touches, filename=filename, ...)
	}
)

setMethod("project", signature(x="SpatRaster"),
	function(x, y, method, mask=FALSE, align=FALSE, gdal=TRUE, res=NULL, origin=NULL, filename="", ...)  {
	
		if (missing(method)) {
			if (is.factor(x)[1] || isTRUE(x@ptr$rgb)) {
				method <- "near"
			} else {
				method <- "bilinear"
			}
		} else {
			method <- method[1]
		}
		if (method == "ngb") {
			method <- "near"
			warn("project", "argument 'method=ngb' is deprecated, it should be 'method=near'")
		}
		opt <- spatOptions(filename, ...)

		if (inherits(y, "SpatRaster")) {
			if (gdal) {
				x@ptr <- x@ptr$warp(y@ptr, "", method, mask[1], align[1], FALSE, opt)
			} else {
				if (align) {
					y <- project(rast(x), y, align=TRUE)
				}
				x@ptr <- x@ptr$resample(y@ptr, method, mask[1], TRUE, opt)
			}
		} else {
			if (!is.character(y)) {
				warn("project,SpatRaster", "argument y (the crs) should be a character value")
				y <- as.character(crs(y))
			}
			if (!is.null(res) || !is.null(origin)) {
				tmp <- project(rast(x), y)
				if (!is.null(res)) res(tmp) <- res
				if (!is.null(origin)) origin(tmp) <- origin
				return(project(x, tmp, method=method, mask=mask, align=align, gdal=gdal, filename=filename, ...))
			}
			if (gdal) {
				x@ptr <- x@ptr$warp(SpatRaster$new(), y, method, mask, FALSE, FALSE, opt)
			} else {
				y <- project(rast(x), y)
				x@ptr <- x@ptr$resample(y@ptr, method, mask[1], TRUE, opt)
			}
		}
		messages(x, "project")
	}
)


setMethod("project", signature(x="SpatVector"),
	function(x, y)  {
		if (!is.character(y)) {
			y <- as.character(crs(y))
		}
		x@ptr <- x@ptr$project(y)
		messages(x, "project")
	}
)

setMethod("project", signature(x="matrix"),
    function(x, from, to)  {
        if (ncol(x) != 2) {
			error("project", "x must have two columns")
		}
        if (missing(from)) {
			error("project", "'from' cannot be missing")
		}
        if (missing(to)) {
			error("project", "'to' cannot be missing")
		}
        if (!is.character(from)) {
           from <- as.character(crs(from))
        }
        if (!is.character(to)) {
           to <- as.character(crs(to))
        }
		v <- vect(x, type="line", crs=from)
        v@ptr <- v@ptr$project(to)
        messages(v, "project")
        crds(v)
    }
)


setMethod("quantile", signature(x="SpatRaster"),
	function(x, probs=seq(0, 1, 0.25), na.rm=FALSE, filename="", ...) {
		opt <- spatOptions(filename, ...)
		x@ptr <- x@ptr$quantile(probs, na.rm[1], opt)
		messages(x, "quantile")
	}
)


setMethod("quantile", signature(x="SpatVector"),
	function(x, probs=seq(0, 1, 0.25), ...) {
		x <- values(x)
		cls <- sapply(x, class)
		i <- cls != "character"
		if (!any(i)) error("quantile", "no numeric variables")
		x <- x[, i, drop=FALSE]
		apply(x, 2, function(i) quantile(i, probs=probs, ...))
	}
)


setMethod("rectify", signature(x="SpatRaster"),
	function(x, method="bilinear", aoi=NULL, snap=TRUE, filename="", ...) {
		opt <- spatOptions(filename, ...)
		if (!is.null(aoi)) {
			if (inherits(aoi, "SpatExtent")) {
				aoi <- rast(aoi)
				useaoi <- 1
			} else if (inherits(aoi, "SpatRaster")) {
				aoi <- rast(aoi)
				useaoi <- 2
			} else {
				error("rectify", "ext must be a SpatExtent or SpatRaster")
			}
		} else {
			aoi <- rast()
			useaoi <- 0
		}
		snap <- as.logical(snap)
		x@ptr <- x@ptr$rectify(method, aoi@ptr, useaoi, snap, opt)
		messages(x, "rectify")
	}
)

setMethod("resample", signature(x="SpatRaster", y="SpatRaster"),
	function(x, y, method, filename="", ...)  {

		if (missing(method)) {
			method <- ifelse(is.factor(x)[1], "near", "bilinear")
		}
		if (method == "ngb") {
			method <- "near"
			warn("resample", "argument 'method=ngb' is deprecated, it should be 'method=near'")
		}
		xcrs = crs(x)
		ycrs = crs(y)
		if ((xcrs == "") && (ycrs != "")) {
			crs(x) <- ycrs
		}
		if ((ycrs == "") && (xcrs != "")) {
			crs(y) <- xcrs
		}
		opt <- spatOptions(filename, ...)
#		if (gdal) {
			x@ptr <- x@ptr$warp(y@ptr, "", method, FALSE, FALSE, TRUE, opt)
#		} else {
#			x@ptr <- x@ptr$resample(y@ptr, method, FALSE, TRUE, opt)
#		}
		messages(x, "resample")
	}
)



setMethod("impose", signature(x="SpatRasterCollection"),
	function(x, y, filename="", ...)  {
		stopifnot(inherits(y, "SpatRaster"))
		opt <- spatOptions(filename, ...)
		r <- rast()
		r@ptr <- x@ptr$morph(y@ptr, opt)
		messages(r, "impose")
	}
)



setMethod("rev", signature(x="SpatRaster"),
	function(x) {
		opt <- spatOptions("", FALSE, list())
		x@ptr <- x@ptr$reverse(opt)
		messages(x, "rev")
	}
)

setMethod("rotate", signature(x="SpatRaster"),
	function(x, left=TRUE, filename="", ...) {
		opt <- spatOptions(filename, ...)
		x@ptr <- x@ptr$rotate(left, opt)
		messages(x, "rotate")
	}
)


setMethod("segregate", signature(x="SpatRaster"),
	function(x, classes=NULL, keep=FALSE, other=0, round=FALSE, digits=0, filename="", ...) {
		opt <- spatOptions(filename, ...)
		if (is.null(classes)) classes <- 1[0]
		x@ptr <- x@ptr$separate(classes, keep[1], other[1], round[1], digits[1], opt)
		messages(x, "segregate")
	}
)


setMethod("shift", signature(x="SpatRaster"),
	function(x, dx=0, dy=0, filename="", ...) {
		opt <- spatOptions(filename, ...)
		x@ptr <- x@ptr$shift(dx, dy, opt)
		messages(x, "shift")
	}
)


setMethod("shift", signature(x="SpatExtent"),
	function(x, dx=0, dy=0) {
		s <- c(dx[1], dx[1], dy[1], dy[1])
		ext(as.vector(x) + s)
	}
)


setMethod("shift", signature(x="SpatVector"),
	function(x, dx=0, dy=0) {
		x@ptr <- x@ptr$shift(dx, dy)
		messages(x, "shift")
	}
)

setMethod("rescale", signature(x="SpatRaster"),
	function(x, fx=0.5, fy=fx, x0, y0) {
		stopifnot(fx > 0)
		stopifnot(fy > 0)
		e <- as.vector(ext(x))
		if (missing(x0)) {
			x0 <- mean(e[1:2])
		}
		if (missing(y0)) {
			y0 <- mean(e[3:4])
		}
		ex = x0 + fx * (e[1:2] - x0);
		ey = y0 + fy * (e[3:4] - y0);
		x@ptr <- x@ptr$deepcopy()
		ext(x) <- ext(c(ex, ey))
		messages(x, "rescale")
	}
)

setMethod("rescale", signature(x="SpatVector"),
	function(x, fx=0.5, fy=fx, x0, y0) {
		stopifnot(fx > 0)
		stopifnot(fy > 0)
		e <- as.vector(ext(x))
		if (missing(x0)) {
			x0 <- mean(e[1:2])
		}
		if (missing(y0)) {
			y0 <- mean(e[3:4])
		}
		x@ptr <- x@ptr$rescale(fx, fy, x0[1], y0[1])
		messages(x, "rescale")
	}
)


setMethod("scale", signature(x="SpatRaster"),
	function(x, center=TRUE, scale=TRUE) {

		opt <- spatOptions()

		if (is.logical(center)) {
			docenter = center[1];
			center = 1[0]
		} else {
			docenter = TRUE
		}
		if (is.logical(scale)) {
			doscale = scale[1]
			scale = 1[0]
		} else {
			doscale = TRUE;
		}
		x@ptr <- x@ptr$scale(center, docenter, scale, doscale, opt)
		messages(x, "scale")
	}
)



setMethod("stretch", signature(x="SpatRaster"),
	function(x, minv=0, maxv=255, minq=0, maxq=1, smin=NA, smax=NA, histeq=FALSE, scale=1, filename="", ...) {
		opt <- spatOptions(filename, ...)
		if (histeq) {
			eqStretch <- function(x){
				ecdfun <- stats::ecdf(x)(x)
				ecdfun(x)
			}
			setValues(x, apply(values(x), 2, function(i) stats::ecdf(i)(i))) * scale
		} else {
			x@ptr <- x@ptr$stretch(minv, maxv, minq, maxq, smin, smax, opt)
			messages(x, "stretch")
		}
	}
)



setMethod("summary", signature(object="SpatRaster"),
	function(object, size=100000, warn=TRUE, ...)  {
		if (!hasValues(object)) {
			warn("summary", "SpatRaster has no values")
			return(invisible())
		}
		if (warn && (ncell(object) > size)) {
			warn("summary", "used a sample")
		}
		s <- spatSample(object, size, method="regular")
		summary(s, ...)
	}
)


setMethod("summary", signature(object="SpatVector"),
	function(object, ...)  {
		summary(as.data.frame(object), ...)
	}
)


setMethod("t", signature(x="SpatRaster"),
	function(x) {
		opt <- spatOptions()
		x@ptr <- x@ptr$transpose(opt)
		messages(x, "t")
	}
)

setMethod("t", signature(x="SpatVector"),
	function(x) {
		x@ptr <- x@ptr$transpose()
		messages(x, "t")
	}
)


setMethod("terrain", signature(x="SpatRaster"),
	function(x, v="slope", neighbors=8, unit="degrees", filename="", ...) {
		unit <- match.arg(unit, c("degrees", "radians"))
		opt <- spatOptions(filename, ...)
		seed <- ifelse("flowdir" %in% v, .seed(), 0)
		x@ptr <- x@ptr$terrain(v, neighbors[1], unit=="degrees", seed, opt)
		messages(x, "terrain")
	}
)


setMethod("trim", signature(x="SpatRaster"),
	function(x, padding=0, value=NA, filename="", ...) {
		opt <- spatOptions(filename, ...)
		x@ptr <- x@ptr$trim(value[1], padding[1], opt)
		messages(x, "trim")
	}
)

setMethod("trans", signature(x="SpatRaster"),
	function(x, filename="", ...) {
		opt <- spatOptions(filename, ...)
		x@ptr <- x@ptr$transpose(opt)
		messages(x, "trans")
	}
)


setMethod("unique", signature(x="SpatRaster", incomparables="ANY"),
	function(x, incomparables=FALSE, na.rm=TRUE, as.raster=FALSE) {

		opt <- spatOptions()

		if (as.raster) incomparables = FALSE
		u <- x@ptr$unique(incomparables, na.rm[1], opt)

		isfact <- is.factor(x)
		if (any(isfact)) {
			ff <- which(isfact)
			levs <- cats(x)
			for (f in ff) {
				lvs <- levs[[f]]
				fv <- factor(lvs[,2])
				i <- match(u[[f]], lvs[,1])
				u[[f]] = fv[i]
			}
		}
		if (!incomparables) {
			#if (!length(u)) return(u)
			u <- do.call(data.frame, u)
			colnames(u) <- names(x)
			if (nrow(u) == 0) {
				if (as.raster) {
					return(deepcopy(x))
				} else {
					return(NULL)
				}
			}
		}

		if (na.rm & (NCOL(u) > 1)) {
			i <- apply(is.na(u), 1, all)
			u <- u[!i, , drop=FALSE]
		}
		if (as.raster) {
			if (any(isfact)) {
				warn("unique", "cannot do 'as.raster=TRUE' with categorical rasters (but you can use 'concats' for that)")
				return(u)
			}
			uid <- 1:nrow(u)
			x <- subst(x, u, uid-1)
			lab <- apply(u, 1, function(i) paste(i, collapse="_"))
			set.cats(x, 1, data.frame(ID=uid-1, label=lab, u), 2)
			return(x)
		}
		u
	}
)


setMethod("unique", signature(x="SpatVector", incomparables="ANY"),
	function(x, incomparables=FALSE, ...) {
		u <- unique(as.data.frame(x, geom="WKT"), incomparables=incomparables, ...)
		vect(u, geom="geometry", crs(x))
	}
)


setMethod("labels", signature(object="SpatRaster"),
	function(object, ...)  {
		names(object)
	}
)


setMethod("sort", signature(x="SpatRaster"),
	function (x, decreasing=FALSE, filename="", ...) {
		opt <- spatOptions(filename, ...)
		x@ptr <- x@ptr$sort(decreasing[1], opt)
		messages(x)
	}
)



#### EC 20210702
setMethod("watershed2", signature(p="SpatRaster",pp_offset="integer"), 
          function(p,pp_offset,filename="", ...) { 
            ##  v <- match.arg(unique(v), c("aspect", "flowdir", "roughness", "slope", "TPI", "TRI"))
            ##  unit <- match.arg(unit, c("degrees", "radians"))
            ##  opt <- spatOptions(filename, ...)
            ##  seed <- ifelse("flowdirection" %in% v, .seed(), 0)
            print("watershed")
            opt <- spatOptions(filename, ...)
            p@ptr <- p@ptr$watershed2(as.integer(pp_offset-1),opt)
            messages(p, "watershed2") ## EC 20210318
            return(p)
            ## p@ptr <- uu
            ##messages(p, "watershed2")
          }
          
)

#### END EC 20210702

#### EC 20220809
setMethod("pitfinder2", signature(p="SpatRaster"), 
          function(p,filename="", ...) { 
      
            opt <- spatOptions(filename, ...)
            p@ptr <- p@ptr$pitfinder2(opt)
            messages(p, "pitfinder2") ## EC 20210318
            return(p)
            ## p@ptr <- uu
            ##messages(p, "watershed2")
          }
          
)
#### END EC 20220809
