# Author: Robert J. Hijmans
# Date : October 2018
# Version 1.0
# License GPL v3

setMethod("is.rotated", signature(x="SpatRaster"),
	function(x) {
		x@pnt$is_rotated()
	}
)


setMethod("rangeFill", signature(x="SpatRaster"),
	function(x, limit, circular=FALSE, filename="", ...) {
		opt <- spatOptions(filename, ...)
		x@pnt <- x@pnt$fill_range(limit, circular, opt)
		messages(x, "rangeFill")
	}
)

setMethod("meta", signature(x="SpatRaster"),
	function(x, layers=FALSE) {
		f <- function(i) {
			if (length(i) == 0) {
				matrix(ncol=2, nrow=0)
			} else {
				matrix(unlist(regmatches(i, regexpr("=", i), invert=TRUE)), ncol=2, byrow=TRUE)
			}
		}
		lapply(x@pnt$metadata(layers), f)
	}
)


setMethod("weighted.mean", signature(x="SpatRaster", w="numeric"),
	function(x, w, na.rm=FALSE, filename="", ...) {
		opt <- spatOptions(filename, ...)
		x@pnt <- x@pnt$wmean_vect(w, na.rm, opt)
		messages(x, "weighted.mean")
	}
)


setMethod("weighted.mean", signature(x="SpatRaster", w="SpatRaster"),
	function(x, w, na.rm=FALSE, filename="", ...) {
		opt <- spatOptions(filename, ...)
		x@pnt <-x@pnt$wmean_rast(w@pnt, na.rm, opt)
		messages(x, "weighted.mean")
	}
)



setMethod("patches", signature(x="SpatRaster"),
	function(x, directions=4, zeroAsNA=FALSE, allowGaps=TRUE, filename="", ...) {
		if (allowGaps) {
			opt <- spatOptions(filename, ...)
			x@pnt <- x@pnt$patches(directions[1], zeroAsNA[1], opt)
			messages(x, "patches")
		} else {
			opt <- spatOptions()
			x@pnt <- x@pnt$patches(directions[1], zeroAsNA[1], opt)
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
		x@pnt$origin
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
		x@pnt <- y@pnt$align(x@pnt, tolower(snap))
		#messages(x, "align")
		x
	}
)

setMethod("align", signature(x="SpatExtent", y="numeric"),
	function(x, y) {
		x@pnt <- x@pnt$align(y, "")
		x
	}
)


setMethod("intersect", signature(x="SpatRaster", y="SpatRaster"),
	function(x, y) {
		opt <- spatOptions() 
		x@pnt <- x@pnt$intersect(y@pnt, opt)
		messages(x)
	}
)


setMethod("cellSize", signature(x="SpatRaster"),
	function(x, mask=FALSE, lyrs=FALSE, unit="m", transform=TRUE, rcx=100, filename="", ...) {
		opt <- spatOptions(filename, ...)
		if (!lyrs) x <- x[[1]]
		x@pnt <- x@pnt$rst_area(mask, unit, transform, rcx, opt)
		messages(x, "cellSize")
	}
)


setMethod("atan2", signature(y="SpatRaster", x="SpatRaster"),
	function(y, x) {
		opt <- spatOptions(filename="", overwrite=TRUE)
		y@pnt <- y@pnt$atan2(x@pnt, opt)
		messages(y, "atan2")
	}
)

setMethod("atan_2", signature(y="SpatRaster", x="SpatRaster"),
	function(y, x, filename="", ...) {
		opt <- spatOptions(filename=filename, ...)
		y@pnt <- y@pnt$atan2(x@pnt, opt)
		messages(y, "atan_2")
	}
)


setMethod("boundaries", signature(x="SpatRaster"),
	function(x, classes=FALSE, inner=TRUE, directions=8, falseval=0, filename="", ...) {
		opt <- spatOptions(filename, ...)
		type <- ifelse(inner[1], "inner", "outer")
		x@pnt <- x@pnt$boundaries(classes[1], type, directions[1], falseval[1], opt)
		messages(x, "boundaries")
	}
)


.collapseSources <- function(x) {
	x@pnt <- x@pnt$collapse_sources()
	messages(x, "tighten")
}

setMethod("deepcopy", signature("SpatRaster"),
	function(x) {
		x@pnt <- x@pnt$deepcopy()
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
		if (x@pnt$same(value@pnt)) {
			x@pnt <- x@pnt$deepcopy()
		}
		opt <- spatOptions()
		x@pnt$addSource(value@pnt, FALSE, opt)
		messages(x, "add")
	}
)

setMethod("tighten", signature("SpatRaster"),
	function(x) {
		x@pnt <- x@pnt$collapse_sources()
		messages(x, "tighten")
	}
)

setMethod("tighten", signature("SpatRasterDataset"),
	function(x) {
		y <- new("SpatRaster")
		y@pnt <- x@pnt$collapse()
		messages(y, "tighten")
	}
)


#setMethod("c", signature(x="SpatRaster"),
#	function(x, ...) {
#		s <- sds(list(x, ...))
#		x@pnt <- s@pnt$collapse()
#		x <- messages(x, "c")
#		try( x@pnt <- x@pnt$collapse_sources() )
#		messages(x, "c")
#	}
#)



#cbind.SpatVector <- function(x, y, ...) {
#	if (inherits(y, "SpatVector")) {
#		y <- y@pnt$df
#	} else {
#		stopifnot(inherits(y, "data.frame"))
#		y <- terra:::.makeSpatDF(y)
#	}
#	x@pnt <- x@pnt$cbind(y)
#	messages(x, "cbind")
#}

cbind.SpatVector <- function(x, y, ...) {
	dots <- list(y, ...)
	for (y in dots) {
		if (inherits(y, "SpatVector")) {
			y <- y@pnt$df
		} else {
			stopifnot(inherits(y, "data.frame"))
			y <- .makeSpatDF(y)
		}
		x@pnt <- x@pnt$cbind(y)
		x <- messages(x, "cbind")
	}
	x
}

rbind.SpatVector <- function(x, y, ...) {
	skipped <- FALSE
	if (missing(y) || is.null(y)) return(x)
	stopifnot(inherits(y, "SpatVector"))
	s <- svc(x, y, ...)
	x@pnt <- s@pnt$append()
	messages(x, "rbind")

#	dots <- list(...)
#	if (is.null(dots)) {
#		x@pnt <- x@pnt$rbind(y@pnt, FALSE)
#	} else {

	#if (!is.null(dots)) {
	#	for (y in dots) {
	#		stopifnot(inherits(y, "SpatVector"))
	#		x@pnt <- x@pnt$rbind(y@pnt, FALSE)
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


setMethod("rep", signature(x="SpatRaster"),
	function(x, ...) {
		i <- rep(1:nlyr(x), ...)
		x[[i]]
	}
)


setMethod("clamp", signature(x="SpatRaster"),
	function(x, lower=-Inf, upper=Inf, values=TRUE, filename="", ...) {
		opt <- spatOptions(filename, ...)
		rlow <- inherits(lower, "SpatRaster")
		rupp <- inherits(upper, "SpatRaster")
		r <- rast()
		if (rlow & rupp) {
			x@pnt <- x@pnt$clamp_raster(lower@pnt, upper@pnt, NA, NA, values[1], opt)
		} else if (rlow) {
			if (any(is.na(upper))) error("clamp", "upper limit cannot be NA")
			x@pnt <- x@pnt$clamp_raster(lower@pnt, r@pnt, NA, upper, values[1], opt)
		} else if (rupp) {
			if (any(is.na(lower))) error("clamp", "lower limit cannot be NA")
			x@pnt <- x@pnt$clamp_raster(r@pnt, upper@pnt, lower, NA, values[1], opt)
		} else {
			if (any(is.na(lower))) error("clamp", "lower limit cannot be NA")
			if (any(is.na(upper))) error("clamp", "upper limit cannot be NA")
			x@pnt <- x@pnt$clamp_raster(r@pnt, r@pnt, lower, upper, values[1], opt)
			#x@pnt <- x@pnt$clamp(lower, upper, values[1], opt)
		}
		messages(x, "clamp")
	}
)



setMethod("clamp_ts", signature(x="SpatRaster"),
	function(x, min=FALSE, max=TRUE, filename="", ...) {
		opt <- spatOptions(filename, ...)
		x@pnt <- x@pnt$clamp_ts(min, max, opt)
		messages(x, "clamp_ts")
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
    x@pnt <- x@pnt$classify(as.vector(rcl), NCOL(rcl), right, include.lowest, others, othersValue, bylayer[1], brackets[1], keepcats, opt)
	messages(x, "classify")
}
)

setMethod("subst", signature(x="SpatRaster"),
function(x, from, to, others=NULL, raw=FALSE, filename="", ...) {
	opt <- spatOptions(filename, ...)

	if (inherits(from, "data.frame")) {
		from <- as.matrix(from)
	}
	if (inherits(to, "data.frame")) {
		to <- as.matrix(to)
	}
	tom <- inherits(to, "matrix")
	frm <- inherits(from, "matrix")
	if (tom && frm) {
		error("subst", "either 'to' or 'from' can be a matrix, not both")
	}

	if (NROW(from) < NROW(to)) {
		error("subst", "from is shorter than to")
	}
	fromc <- inherits(from[1], "character")
	toc <- inherits(to[1], "character")
	if (raw && fromc) {
		error("subst", "if 'raw=TRUE', 'from' cannot have character values")
	}
	keepcats <- FALSE
	if (any(is.factor(x))) {
		if (nlyr(x) > 1) {
			error("subst", "you can only use 'subst' with categorical layers if x has a single layer")
		}
		if (inherits(to, "matrix")) {
			if (ncol(to) > 1) {
				warn("subst", "only the first column of 'to' is used with factors")
			}
			to <- as.vector(to[,1])
		}
		levs <- levels(x)[[1]]
		if (!raw) {
			from <- levs[,1][match(from, levs[,2])]
			if (all(is.na(from))) {
				warn("subst", "all 'from' values are missing, returning a copy of 'x'")
				return(deepcopy(x))
			}
		}
		i <- is.na(from)
		if (any(i)) {
			to <- rep_len(to, length(from))
			from <- from[!i]
			to <- to[!i]
		}
		if (!raw) {
			toto <- levs[,1][match(to, levs[,2])]
			if (any(is.na(toto))) { # add new levels
				i <- which(is.na(toto))
				m <- cbind(max(levs[,1]) + 1:length(i), to[i])
				colnames(m) <- colnames(levs)
				levs <- rbind(levs, m)
				to <- levs[,1][match(to, levs[,2])]
				levels(x) <- levs
			} else {
				to <- toto
			}
		}
		keepcats <- TRUE
	} else if (fromc || toc) {
		error("subst", "from or to has character values but x is not categorical")
	}
	
	if (is.null(others)) {
		setothers <- FALSE
		others <- NA
	} else {
		setothers <- TRUE
		others <- others[1]
	}
	
	if (tom) {
		nms <- colnames(to)
		if (!is.null(nms)) 	opt$names = nms
		x@pnt <- x@pnt$replaceValues(from, to, ncol(to), setothers, others, keepcats, opt)
	} else if (frm) {
		x@pnt <- x@pnt$replaceValues(as.vector(t(from)), to, -ncol(from), setothers, others, FALSE, opt)
	} else {
		x@pnt <- x@pnt$replaceValues(from, to, 0, setothers, others, keepcats, opt)
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
	function(x, y, snap="near", mask=FALSE, touches=TRUE, extend=FALSE, filename="", ...) {
		opt <- spatOptions(filename, ...)
		if (mask) {
			if (inherits(y, "SpatVector")) {
				x@pnt <- x@pnt$crop_mask(y@pnt, snap[1], touches[1], extend[1], opt)
			} else if (inherits(y, "sf")) {
				y <- vect(y)
				x@pnt <- x@pnt$crop_mask(y@pnt, snap[1], touches[1], extend[1], opt)
			} else if (inherits(y, "SpatRaster")) {
				mopt <- spatOptions(filename="", ...)
				e <- ext(y)
				x@pnt <- x@pnt$crop(e@pnt, snap[1], extend[1], mopt)
				x <- messages(x, "crop")
				if (compareGeom(x, y, lyrs=FALSE, crs=FALSE, warncrs=FALSE, ext=TRUE, rowcol=TRUE, res=FALSE, stopOnError=FALSE, messages=FALSE)) {
					return(mask(x, y, filename=filename, ...))
				} else {
					warn("crop", "cannot mask with a raster that is not aligned to x")
					# should check earlier if masking will be possible to avoid writing to disk twice.
					if (filename != "") {
						x <- writeRaster(x, filename, ...)
					}
					return(x)
				}
			} else {
				y <- .getExt(y, method="crop")
				warn("crop", paste("mask=TRUE is ignored for", class(y)[1], "objects"))
				x@pnt <- x@pnt$crop(y@pnt, snap[1], extend[1], opt)
			}
		} else {
			y <- .getExt(y, method="crop")
			x@pnt <- x@pnt$crop(y@pnt, snap[1], extend[1], opt)
		}
		messages(x, "crop")
	}
)


setMethod("crop", signature(x="SpatRasterDataset", y="ANY"),
	function(x, y, snap="near", extend=FALSE) {
		opt <- spatOptions()
		y <- .getExt(y, method="crop")
		x@pnt <- x@pnt$crop(y@pnt, snap[1], extend[1], opt)
		messages(x, "crop")
	}
)

setMethod("crop", signature(x="SpatRasterCollection", y="ANY"),
	function(x, y, snap="near", extend=FALSE) {
		y <- .getExt(y, method="crop")
		opt <- spatOptions()
		x@pnt <- x@pnt$crop(y@pnt, snap[1], extend[1], double(), opt)
		messages(x, "crop")
	}
)


setMethod("selectRange", signature(x="SpatRaster"),
	function(x, y, z=1, repint=0, filename="", ...) {
		opt <- spatOptions(filename, ...)
		x@pnt <- x@pnt$selRange(y@pnt, z, repint, opt)
		messages(x, "selectRange")
	}
)

setMethod("cover", signature(x="SpatRaster", y="SpatRaster"),
	function(x, y, values=NA, filename="", ...) {
		opt <- spatOptions(filename, ...)
		x@pnt <- x@pnt$cover(y@pnt, values, opt)
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
		x@pnt <- x@pnt$arith_rast(y@pnt, "-", FALSE, opt)
		messages(x, "diff")
	}
)


setMethod("disagg", signature(x="SpatRaster"),
	function(x, fact, method="near", filename="", ...) {
		method <- match.arg(tolower(method), c("near", "bilinear"))
		if (method == "bilinear") {
			y <- disagg(rast(x), fact)
			r <- resample(x, y, "bilinear", filename=filename, ...)
			return(r)
		}
		opt <- spatOptions(filename, ...)
		x@pnt <- x@pnt$disaggregate(fact, opt)
		messages(x, "disagg")
	}
)


setMethod("flip", signature(x="SpatRaster"),
	function(x, direction="vertical", filename="", ...) {
		d <- match.arg(direction, c("vertical", "horizontal"))
		opt <- spatOptions(filename, ...)
		x@pnt <- x@pnt$flip(d == "vertical", opt)
		messages(x, "flip")
	}
)



setMethod("mask", signature(x="SpatRaster", mask="SpatRaster"),
	function(x, mask, inverse=FALSE, maskvalues=NA, updatevalue=NA, filename="", ...) {
		opt <- spatOptions(filename, ...)
		x@pnt <- x@pnt$mask_raster(mask@pnt, inverse[1], maskvalues, updatevalue[1], opt)
		messages(x, "mask")
	}
)

setMethod("mask", signature(x="SpatRaster", mask="SpatVector"),
	function(x, mask, inverse=FALSE, updatevalue=NA, touches=TRUE, filename="", ...) {
		opt <- spatOptions(filename, ...)
		x@pnt <- x@pnt$mask_vector(mask@pnt, inverse[1], updatevalue[1], touches[1], opt)
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
	function(x, y, method, mask=FALSE, align=FALSE, gdal=TRUE, res=NULL, origin=NULL, threads=FALSE, filename="", ..., by_util = FALSE)  {

		if (missing(method)) {
			if (is.factor(x)[1] || isTRUE(x@pnt$rgb)) {
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
		opt <- spatOptions(filename, threads=threads, ...)

		if (inherits(y, "SpatRaster")) {
			if (gdal) {

				if (by_util) {
						x@pnt <- x@pnt$warp_by_util(y@pnt, "", method, mask[1], align[1], FALSE, opt)
				} else {
					x@pnt <- x@pnt$warp(y@pnt, "", method, mask[1], align[1], FALSE, opt)
				}
			} else {
				if (align) {
					y <- project(rast(x), y, align=TRUE)
				}
				x@pnt <- x@pnt$resample(y@pnt, method, mask[1], TRUE, opt)
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

				if (by_util) {
					x@pnt <- x@pnt$warp_by_util(SpatRaster$new(), y, method, mask, FALSE, FALSE, opt)
					
				} else {
					x@pnt <- x@pnt$warp(SpatRaster$new(), y, method, mask, FALSE, FALSE, opt)
				}
			} else {
				y <- project(rast(x), y)
				x@pnt <- x@pnt$resample(y@pnt, method, mask[1], TRUE, opt)
			}
		}
		messages(x, "project")
	}
)


setMethod("project", signature(x="SpatVector"),
	function(x, y, partial=FALSE)  {
		if (!is.character(y)) {
			y <- as.character(crs(y))
		}
		x@pnt <- x@pnt$project(y, partial)
		messages(x, "project")
	}
)

setMethod("project", signature(x="SpatExtent"),
	function(x, from, to)  {
		if (missing(from)) error("project", "'from' cannot be missing")
		if (missing(to)) error("project", "'to' cannot be missing")
		x <- as.polygons(x, crs=from)
		x <- densify(x, 10000)
		ext(project(x, to))
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
		#v <- vect(x, type="line", crs=from)
        #v@pnt <- v@pnt$project(to)
        v <- vect()
		xy <- v@pnt$project_xy(x[,1], x[,2], from, to)
		messages(v, "project")
		matrix(xy, ncol=2)
    }
)


setMethod("quantile", signature(x="SpatRaster"),
	function(x, probs=seq(0, 1, 0.25), na.rm=FALSE, filename="", ...) {
		opt <- spatOptions(filename, ...)
		x@pnt <- x@pnt$quantile(probs, na.rm[1], opt)
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
		x@pnt <- x@pnt$rectify(method, aoi@pnt, useaoi, snap, opt)
		messages(x, "rectify")
	}
)

setMethod("resample", signature(x="SpatRaster", y="SpatRaster"),
	function(x, y, method, threads=FALSE, filename="", ...)  {

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
		opt <- spatOptions(filename, threads=threads, ...)
#		if (gdal) {
			x@pnt <- x@pnt$warp(y@pnt, "", method, FALSE, FALSE, TRUE, opt)
#		} else {
#			x@pnt <- x@pnt$resample(y@pnt, method, FALSE, TRUE, opt)
#		}
		messages(x, "resample")
	}
)



setMethod("impose", signature(x="SpatRasterCollection"),
	function(x, y, filename="", ...)  {
		stopifnot(inherits(y, "SpatRaster"))
		opt <- spatOptions(filename, ...)
		r <- rast()
		r@pnt <- x@pnt$morph(y@pnt, opt)
		messages(r, "impose")
	}
)



setMethod("rev", signature(x="SpatRaster"),
	function(x) {
		opt <- spatOptions("", FALSE, list())
		x@pnt <- x@pnt$reverse(opt)
		messages(x, "rev")
	}
)

setMethod("rotate", signature(x="SpatRaster"),
	function(x, left=TRUE, filename="", ...) {
		opt <- spatOptions(filename, ...)
		x@pnt <- x@pnt$rotate(left, opt)
		messages(x, "rotate")
	}
)

setMethod("rotate", signature(x="SpatVector"),
	function(x, longitude=0, split=FALSE, left=TRUE, normalize=FALSE) {
		if (split) {
			e <- ext(x)
			if ((longitude < e$xmin) || (longitude > e$xmax)) {
				if (left) {
					return(shift(x, -360))
				} else {
					return(shift(x, 360))
				}
			}
			e <- as.vector(floor(e) + 1)
			ew <- ext(c(e[1], longitude, e[3:4]))
			ee <- ext(c(longitude, e[2:4]))
			x$unique_id_for_aggregation <- 1:nrow(x)
			xw <- crop(x, ew)
			xe <- crop(x, ee)
			if (left) {
				xe <- shift(xe, -360)
			} else {
				xw <- shift(xw, 360)
			}
			out <- rbind(xe, xw)
			if (nrow(out) > nrow(x)) {
				out <- aggregate(out, "unique_id_for_aggregation", id=F)
				i <- match(out$unique_id_for_aggregation, x$unique_id_for_aggregation)
				values(out) <- values(x)[i,,drop=FALSE]
				x <- out
			} 
			x$unique_id_for_aggregation <- NULL
		} else {
			x@pnt <- x@pnt$rotate_longitude(longitude, left)
			x <- messages(x, "rotate")
		}
		if (normalize) {
			x <- normalize.longitude(x)
		}
		x
	}
)

setMethod("segregate", signature(x="SpatRaster"),
	function(x, classes=NULL, keep=FALSE, other=0, round=FALSE, digits=0, filename="", ...) {
		opt <- spatOptions(filename, ...)
		if (is.null(classes)) classes <- 1[0]
		x@pnt <- x@pnt$separate(classes, keep[1], other[1], round[1], digits[1], opt)
		messages(x, "segregate")
	}
)


setMethod("shift", signature(x="SpatRaster"),
	function(x, dx=0, dy=0, filename="", ...) {
		opt <- spatOptions(filename, ...)
		x@pnt <- x@pnt$shift(dx, dy, opt)
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
		x@pnt <- x@pnt$shift(dx, dy)
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
		x@pnt <- x@pnt$deepcopy()
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
		x@pnt <- x@pnt$rescale(fx, fy, x0[1], y0[1])
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
		x@pnt <- x@pnt$scale(center, docenter, scale, doscale, opt)
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
			x@pnt <- x@pnt$stretch(minv, maxv, minq, maxq, smin, smax, opt)
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
		s <- spatSample(object, size, method="regular", warn=FALSE)
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
		x@pnt <- x@pnt$transpose(opt)
		messages(x, "t")
	}
)

setMethod("t", signature(x="SpatVector"),
	function(x) {
		x@pnt <- x@pnt$transpose()
		messages(x, "t")
	}
)


setMethod("terrain", signature(x="SpatRaster"),
	function(x, v="slope", neighbors=8, unit="degrees", filename="", ...) {
		unit <- match.arg(unit, c("degrees", "radians"))
		opt <- spatOptions(filename, ...)
		seed <- ifelse("flowdir" %in% v, .seed(), 0)
		x@pnt <- x@pnt$terrain(v, neighbors[1], unit=="degrees", seed, opt)
		messages(x, "terrain")
	}
)


setMethod("viewshed", signature(x="SpatRaster"),
	function(x, loc, observer=1.80, target=0, curvcoef=6/7, output="yes/no", filename="", ...) {
		opt <- spatOptions(filename, ...)
		z <- rast()
		if (length(loc) == 1) {
			loc <- xyFromCell(x, loc)
		}
		outops <- c("yes/no", "sea", "land")
		output <- match.arg(tolower(output), outops)
		output <- match(output, outops)
		z@pnt <- x@pnt$view(c(loc[1:2], observer[1], target[1]), c(1,0,2,3), curvcoef, 2, 0, output, opt)
		messages(z, "viewshed")
	}
)

setMethod("sieve", signature(x="SpatRaster"),
	function(x, threshold, directions=8, filename="", ...) {
		opt <- spatOptions(filename, ...)
		x@pnt <- x@pnt$sieve(threshold[1], directions[1], opt)
		messages(x, "sieve")
	}
)


setMethod("trim", signature(x="SpatRaster"),
	function(x, padding=0, value=NA, filename="", ...) {
		opt <- spatOptions(filename, ...)
		padding <- round(padding[1])
		if (padding < 0) {
			error("trim", "padding must be a non-negative integer")
		}
		x@pnt <- x@pnt$trim(value[1], padding, opt)
		messages(x, "trim")
	}
)

setMethod("trans", signature(x="SpatRaster"),
	function(x, filename="", ...) {
		opt <- spatOptions(filename, ...)
		x@pnt <- x@pnt$transpose(opt)
		messages(x, "trans")
	}
)


setMethod("unique", signature(x="SpatRaster", incomparables="ANY"),
	function(x, incomparables=FALSE, digits=NA, na.rm=TRUE, as.raster=FALSE) {

		opt <- spatOptions()

		if (as.raster) incomparables = FALSE
		u <- x@pnt$unique(incomparables, digits, na.rm[1], opt)

		if (!as.raster) {
			u <- get_labels(x, u)
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

		if ((!incomparables) && (na.rm || as.raster)) {
			i <- rowSums(is.na(u)) < ncol(u)
			u <- u[i, , drop=FALSE]
		}

		if (as.raster) {
			lab <- apply(get_labels(x, u), 1, function(i) paste(i, collapse="_"))
			if (!is.na(digits)) {
				x <- round(x, digits)
			} else {
				levels(x) <- NULL
			}
			uid <- 1:nrow(u)		
			x <- subst(x, u, uid-1)
			set.cats(x, 1, data.frame(ID=uid-1, label=lab, u))
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


setMethod("scoff", signature(x="SpatRaster"),
	function(x) {
		out <- x@pnt$getScaleOffset()
		names(out) <- c("scale", "offset")
		do.call(cbind, out)
	}
)

setMethod("scoff<-", signature("SpatRaster"),
	function(x, value) {
		if (is.null(value)) {
			x@pnt <- x@pnt$deepcopy()
			x@pnt$setScaleOffset(1, 0)
		} else {
			if (NCOL(value) != 2) {
				error("scoff<-", "value must be a 2-column matrix")
			}
			x@pnt <- x@pnt$deepcopy()
			value[is.na(value[,1]),1] <- 1
			value[is.na(value[,2]),2] <- 0
			x@pnt$setScaleOffset(value[,1], value[,2])
			x@pnt$setValueType(0)
		}
		messages(x, "scoff<-")
	}
)

setMethod("sort", signature(x="SpatRaster"),
	function (x, decreasing=FALSE, order=FALSE, filename="", ...) {
		opt <- spatOptions(filename, ...)
		x@pnt <- x@pnt$sort(decreasing[1], order[1], opt)
		messages(x, "sort")
	}
)


setMethod("sort", signature(x="SpatVector"),
	function (x, v, decreasing=FALSE) {
		if (is.logical(v)) {
			tmp <- v
			v <- decreasing 
			decreasing <- tmp
		}
		if (length(v) > 1) {
			v <- data.frame(x)[,v]
			i <- do.call(order, lapply(v, function(i) i))
		} else {
			i <- order(x[[v]][[1]], decreasing=decreasing)
		}
		x[i, ]
	}
)


setMethod("sort", signature(x="data.frame"),
	function (x, v, decreasing=FALSE) {
		if (length(v) > 1) {
			v <- data.frame(x)[, v]
			i <- do.call(order, lapply(v, function(i) i))
		} else {
			i <- order(x[[v]], decreasing=decreasing)
		}
		x[i, ]
	}
)
