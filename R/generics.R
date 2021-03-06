# Author: Robert J. Hijmans
# Date : October 2018
# Version 1.0
# License GPL v3


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
	function(x, directions=4, zeroAsNA=FALSE, filename="", ...) {
		opt <- spatOptions(filename, ...)
		x@ptr <- x@ptr$patches(directions[1], zeroAsNA[1], opt)
		messages(x, "patches")
	}
)


setMethod("origin", signature(x="SpatRaster"), 
	function(x) {
		x@ptr$origin
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
	function(x, mask=TRUE, unit="m", transform=TRUE, filename="", ...) {
		opt <- spatOptions(filename, ...)
		x@ptr <- x@ptr$rst_area(mask, unit, transform, opt)
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

setMethod("copy", signature("SpatRaster"), 
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
		x@ptr$addSource(value@ptr)
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
	stopifnot(inherits(y, "SpatVector"))
	x@ptr <- x@ptr$rbind(y@ptr, FALSE)
	x <- messages(x, "rbind")
	dots <- list(...)
	if (!is.null(dots)) {
		for (y in dots) {
			stopifnot(inherits(y, "SpatVector"))
			x@ptr <- x@ptr$rbind(y@ptr, FALSE)
			x <- messages(x, "rbind")
		}
	}
	x
}


setMethod("c", signature(x="SpatRaster"), 
	function(x, ...) {
		skipped <- FALSE
		dots <- list(...)
		for (i in dots) {
			if (inherits(i, "SpatRaster")) {
				x@ptr <- x@ptr$combineSources(i@ptr)
				if (x@ptr$messages$has_error) {
					messages(x, "c")
					return()
				}
			} else {
				skipped = TRUE
			}
		}
		if (skipped) warn("c,SpatRaster", "skipped object that are not SpatRaster")
		messages(x, "c")
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
		x@ptr <- x@ptr$clamp(lower, upper, values[1], opt)
		messages(x, "clamp")
	}
)


setMethod("classify", signature(x="SpatRaster"), 
function(x, rcl, include.lowest=FALSE, right=TRUE, othersNA=FALSE, filename="", ...) {

	bylayer = FALSE
	
	if (is.data.frame(rcl)) {
		rcl <- as.matrix(rcl)
	}

	right <- ifelse(is.na(right), 2, ifelse(right, 1, 0))
	include.lowest <- as.logical(include.lowest[1])

	opt <- spatOptions(filename, ...)
    x@ptr <- x@ptr$classify(as.vector(rcl), NCOL(rcl), right[1], include.lowest[1], othersNA[1], bylayer[1], opt)
	messages(x, "classify")
}
)

setMethod("subst", signature(x="SpatRaster"), 
function(x, from, to, filename="", ...) {
	opt <- spatOptions(filename, ...)
	if (inherits(to, "data.frame")) {
		to <- as.matrix(to)
	}
	if (inherits(to, "matrix")) {
		opt$names = colnames(to)
		x@ptr <- x@ptr$replaceValues(from, to, ncol(to), opt)	
	} else {
		x@ptr <- x@ptr$replaceValues(from, to, -1, opt)
	}
	messages(x, "subst")
}
)


.getExt <- function(x) {
	return(x)
}

setMethod("crop", signature(x="SpatRaster", y="ANY"), 
	function(x, y, snap="near", filename="", ...) {
		opt <- spatOptions(filename, ...)

		if (!inherits(y, "SpatExtent")) {
			e <- try(ext(y), silent=TRUE)
			if (class(e) == "try-error") { 
				e <- try(ext(vect(y)), silent=TRUE)
				if (class(e) == "try-error") { 
					error("crop", "cannot get a SpatExtent from y")
				}
			}
			y <- e
		}
		
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


setMethod("disaggregate", signature(x="SpatRaster"), 
	function(x, fact, method="near", filename="", ...) {
		stopifnot(method %in% c("near", "bilinear"))
		if (method == "bilinear") {
			y <- disaggregate(rast(x), fact)
			r <- resample(x, y, "bilinear", filename=filename, ...)
			return(r)
		}
		opt <- spatOptions(filename, ...)
		x@ptr <- x@ptr$disaggregate(fact, opt)
		messages(x, "disaggregate")
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
	function(x, digits=0, value=NULL, bylayer=TRUE) {

		opt <- spatOptions()

		if (!is.null(value)) {
			if (is.na(digits)) {
				v <- x@ptr$count(value, bylayer[1], FALSE, 0, opt)
			} else {
				v <- x@ptr$count(value, bylayer[1], TRUE, digits, opt)
				value <- round(value, digits)
			}
			if (bylayer) {
				v <- cbind(layer=1:nlyr(x), value=value, count=v)
			} else {
				v <- cbind(value=value, count=v)
			}
		} else {
			if (is.na(digits)) {
				v <- x@ptr$freq(bylayer[1], FALSE, 0, opt)
			} else {
				v <- x@ptr$freq(bylayer[1], TRUE, digits, opt)
			}
			if (bylayer) {
				v <- lapply(1:length(v), function(i) cbind(i, matrix(v[[i]], ncol=2)))
				v <- do.call(rbind, v)
				colnames(v) <- c("layer", "value", "count")
			} else {
				v <- matrix(v[[1]], ncol=2, dimnames=list(NULL, c("value", "count")))
			}
			ff <- is.factor(x)
			if (any(ff)) {
				if (bylayer) {
					v <- data.frame(v)
					v$label <- ""
					f <- which(ff)
					levs <- levels(x)
					for (i in f) {
						g <- levs[[i]]
						k <- v$layer==i
						v$label[k] <- g[v$value[k] + 1]
					}
				} else if (nlyr(x) == 1) {				
					v <- data.frame(v)
					g <- levels(x)[[1]]
					v$label <- g[v$value + 1]
				}
			}
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


setMethod("project", signature(x="SpatRaster"), 
	function(x, y, method="bilinear", mask=FALSE, filename="", ...)  {

		method <- ifelse(method == "ngb", "near", method)
		opt <- spatOptions(filename, ...)
		if (inherits(y, "SpatRaster")) {
			#x@ptr <- x@ptr$warp(y@ptr, method, opt)
			x@ptr <- x@ptr$warp(y@ptr, "", method, mask, opt)
		} else {
			if (!is.character(y)) {
				warn("project,SpatRaster", "crs should be a character value")
				y <- as.character(crs(y))
			}
			#x@ptr <- x@ptr$warpcrs(y, method, opt)
			x@ptr <- x@ptr$warp(SpatRaster$new(), y, method, mask, opt)
		}
		messages(x, "project")
	}
)


setMethod("project", signature(x="SpatVector"), 
	function(x, y)  {
		if (!is.character(y)) {
			y <- crs(y)
		}
		x@ptr <- x@ptr$project(y)
		messages(x, "project")
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
	function(x, y, method="bilinear", filename="", ...)  {
		method <- ifelse(method == "ngb", "near", method)
		opt <- spatOptions(filename, ...)
		x@ptr <- x@ptr$warp(y@ptr, "", method, FALSE, opt)
		messages(x, "resample")
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
	function(x, classes=NULL, keep=FALSE, other=0, filename="", ...) {
		opt <- spatOptions(filename, ...)
		if (is.null(classes)) classes <- 1[0]
		x@ptr <- x@ptr$separate(classes, keep, other, opt)
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
	function(x, minv=0, maxv=255, minq=0, maxq=1, smin=NA, smax=NA, filename="", ...) {
		opt <- spatOptions(filename, ...)
		x@ptr <- x@ptr$stretch(minv, maxv, minq, maxq, smin, smax, opt)
		messages(x, "stretch")
	}
)



setMethod("summary", signature(object="SpatRaster"), 
	function(object, size=100000, warn=TRUE, ...)  {
		if (warn && (ncell(object) > size)) {
			warn("summary", "used a sample")
		}
		summary(spatSample(object, size, method="regular"), ...)
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
		v <- match.arg(unique(v), c("aspect", "flowdir", "roughness", "slope", "TPI", "TRI"))
		unit <- match.arg(unit, c("degrees", "radians"))
		opt <- spatOptions(filename, ...)
		seed <- ifelse("flowdirection" %in% v, .seed(), 0)
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
	function(x, incomparables=FALSE) {
		opt <- spatOptions()
		u <- x@ptr$unique(incomparables, opt)
		if (!incomparables) {
			if (!length(u)) return(u)
			u <- do.call(cbind, u)
			colnames(u) = names(x)
		}
		u
	}
)

setMethod("unique", signature(x="SpatVector", incomparables="ANY"), 
	function(x, incomparables=FALSE, ...) {
		u <- unique(as.data.frame(x, geom="WKT"), incomparables=incomparables, ...)
		v <- vect(u, geom="geometry")
		v$geometry <- NULL
		crs(v) <- crs(x)
		v
	}
)


#setMethod("warp", signature(x="SpatRaster", y="SpatRaster"), 
#	function(x, y, method="bilinear", filename="", overwrite=FALSE, ...)  {
#		opt <- spatOptions(filename, overwrite, ...)
#		x@ptr <- x@ptr$warp(y@ptr, "", method, opt)
#		messages(x, "warp")
#	}
#)

