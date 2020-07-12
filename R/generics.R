# Author: Robert J. Hijmans
# Date : October 2018
# Version 1.0
# License GPL v3

setMethod("length", signature(x="SpatRaster"), 
	function(x) {
		ncell(x)
	}
)	

setMethod("origin", signature(x="SpatRaster"), 
	function(x, ...) {
		x@ptr$origin
	}
)	

setMethod("rectify", signature(x="SpatRaster"), 
	function(x, method="bilinear", filename="", overwrite=FALSE, wopt=list(), ...) {
		opt <- .runOptions(filename, overwrite, wopt)
		x@ptr <- x@ptr$rectify(method, opt)
		show_messages(x, "rectify")
	}
)

setMethod("adjacent", signature(x="SpatRaster"), 
	function(x, cells, directions="rook", include=FALSE, ...) {
		v <- x@ptr$adjacent(cells-1, directions, include)
		show_messages(x, "adjacent")
		v <- do.call(rbind, v)
		return(v+1)
	}
)


setMethod("align", signature(x="SpatExtent", y="SpatRaster"), 
	function(x, y, snap="near", ...) {
		x@ptr <- y@ptr$align(x@ptr, tolower(snap))
		#show_messages(x, "align")
		x
	}
)


setMethod("area", signature(x="SpatRaster"), 
	function(x, sum=TRUE, filename="", overwrite=FALSE, wopt=list(), ...) {
		if (sum) {
			byvalue = FALSE
			if (byvalue) {
				v <- x@ptr$area_by_value()
				v <- lapply(1:length(v), function(i) cbind(i, matrix(v[[i]], ncol=2)))
				v <- do.call(rbind, v)
				colnames(v) <- c("layer", "value", "area")
				return(v)
			} else {
				opt <- .getOptions()
				x@ptr$sum_area(opt)
			}
		} else {
			opt <- .runOptions(filename, overwrite, wopt)
			x@ptr <- x@ptr$rst_area(opt)
			show_messages(x, "area")
		} 
	}
)



setMethod("atan2", signature(y="SpatRaster", x="SpatRaster"),
	function(y, x) { 
		opt <- .runOptions(filename="", overwrite=TRUE, wopt=list())
		y@ptr <- y@ptr$atan2(x@ptr, opt)
		show_messages(y, "atan2")
	}
)
	

setMethod("boundaries", signature(x="SpatRaster"), 
	function(x, classes=FALSE, type="inner", directions=8, filename="", overwrite=FALSE, wopt=list(), ...) {
		opt <- .runOptions(filename, overwrite,wopt)
		x@ptr <- x@ptr$boundaries(classes[1], type[1], directions[1], opt)
		show_messages(x, "boundaries")
	}
)


.collapse <- function(x) {
	x@ptr <- x@ptr$collapse_sources()
	show_messages(x, "collapse")
}

setMethod("c", signature(x="SpatRaster"), 
	function(x, ...) {
		s <- sds(list(x, ...))
		x@ptr <- s@ptr$collapse()
		x <- show_messages(x, "c")		
		x@ptr <- x@ptr$collapse_sources()
		show_messages(x, "c")		
	}
)

#setMethod("c", signature(x="SpatRaster"), 
#	function(x, ...) {
#		dots <- list(...)
#		for (i in dots) {
#			if (inherits(i, "SpatRaster")) {
#				x@ptr <- x@ptr$combineSources(i@ptr)
#			}
#		}
#		show_messages(x, "c")		
#	}
#)


setMethod("rep", signature(x="SpatRaster"), 
	function(x, ...) {
		i <- rep(1:nlyr(x), ...)
		x[[i]]
	}
)


setMethod("clamp", signature(x="SpatRaster"), 
	function(x, lower=-Inf, upper=Inf, values=TRUE, filename="", overwrite=FALSE, wopt=list(), ...) {
		opt <- .runOptions(filename, overwrite,wopt)
		x@ptr <- x@ptr$clamp(lower, upper, values[1], opt)
		show_messages(x, "clamp")
	}
)


setMethod("classify", signature(x="SpatRaster"), 
function(x, rcl, include.lowest=FALSE, right=TRUE, othersNA=FALSE, filename="", overwrite=FALSE, wopt=list(), ...) {
	
	if (is.data.frame(rcl)) {
		rcl <- as.matrix(rcl)
	}

	right <- ifelse(is.na(right), 2, ifelse(right, 1, 0))
	include.lowest <- as.logical(include.lowest[1])

	opt <- .runOptions(filename, overwrite, wopt)
    x@ptr <- x@ptr$classify(as.vector(rcl), NCOL(rcl), right, include.lowest, othersNA, opt)
	show_messages(x, "classify")
}
)


.getExt <- function(x) {
	return(x)
}

setMethod("crop", signature(x="SpatRaster", y="ANY"), 
	function(x, y, snap="near", filename="", overwrite=FALSE, wopt=list(), ...) {
		opt <- .runOptions(filename, overwrite, wopt)

		if (!inherits(y, "SpatExtent")) {
			e <- try(ext(y), silent=TRUE)
			if (class(e) == "try-error") { 
				e <- try(raster::extent(y), silent=TRUE)
				if (class(e) == "try-error") { 
					stop("cannot get an extent from y")
				}
				e <- ext(as.vector(t(raster::bbox(e))))
			}
			y <- e
		}

		x@ptr <- x@ptr$crop(y@ptr, snap[1], opt)
		show_messages(x, "crop")		
	}
)



setMethod("selectRange", signature(x="SpatRaster"), 
	function(x, y, z=1, repint=0, filename="", overwrite=FALSE, wopt=list(), ...) { 
		opt <- .runOptions(filename, overwrite, wopt)
		x@ptr <- x@ptr$selRange(y@ptr, z, repint, opt)
		show_messages(x, "selectRange")		
	}
)

setMethod("cover", signature(x="SpatRaster", y="SpatRaster"), 
	function(x, y, value=NA, filename="", overwrite=FALSE, wopt=list(), ...) {
		opt <- .runOptions(filename, overwrite, wopt)
		x@ptr <- x@ptr$cover(y@ptr, value[1], opt)
		show_messages(x, "cover")		
	}
)


setMethod("diff", signature(x="SpatRaster"), 
	function(x, filename="", overwrite=FALSE, wopt=list(), ...) { 
		n = nlyr(x)
		if (n<2) return(rast(x))
		y = x[[-1]]
		x = x[[-n]]
		opt <- .runOptions(filename, overwrite, wopt)
		x@ptr <- x@ptr$arith_rast(y@ptr, "-", opt)
		show_messages(x, "diff")
	}
)


setMethod("disaggregate", signature(x="SpatRaster"), 
	function(x, fact, method="near", filename="", overwrite=FALSE, wopt=list(), ...) {
		stopifnot(method %in% c("near", "bilinear"))
		if (method == "bilinear") {
			y <- disaggregate(rast(x), fact)
			r <- resample(x, y, "bilinear", filename=filename, overwrite=overwrite, wopt=wopt, ...)
			return(r)
		}
		opt <- .runOptions(filename, overwrite, wopt)
		x@ptr <- x@ptr$disaggregate(fact, opt)
		show_messages(x, "disaggregate")
	}
)


setMethod("flip", signature(x="SpatRaster"), 
	function(x, vertical=TRUE, filename="", overwrite=FALSE, wopt=list(), ...) {
		opt <- .runOptions(filename, overwrite, wopt)
		x@ptr <- x@ptr$flip(vertical, opt)
		show_messages(x, "flip")
	}
)


setMethod("freq", signature(x="SpatRaster"), 
	function(x, bylayer=TRUE, ...) {
		opt <- .runOptions("", TRUE, list())
		v <- x@ptr$freq(bylayer[1], opt)
		if (bylayer) {
			v <- lapply(1:length(v), function(i) cbind(i, matrix(v[[i]], ncol=2)))
			v <- do.call(rbind, v)
			colnames(v) <- c("layer", "value", "count")
		} else {
			v <- matrix(v[[1]], ncol=2, dimnames=list(NULL, c("value", "count")))
		}
		v
	}
)


setMethod("mask", signature(x="SpatRaster", mask="SpatRaster"), 
	function(x, mask, inverse=FALSE, maskvalue=NA, updatevalue=NA, filename="", overwrite=FALSE, wopt=list(), ...) { 
		opt <- .runOptions(filename, overwrite,wopt)
		x@ptr <- x@ptr$mask_raster(mask@ptr, inverse[1], maskvalue[1], updatevalue[1], opt)
		show_messages(x, "mask")		
	}
)

setMethod("mask", signature(x="SpatRaster", mask="SpatVector"), 
	function(x, mask, inverse=FALSE, updatevalue=NA, touches=is.lines(mask), filename="", overwrite=FALSE, wopt=list(), ...) { 
		opt <- .runOptions(filename, overwrite,wopt)
		x@ptr <- x@ptr$mask_vector(mask@ptr, inverse[1], updatevalue[1], opt)
		show_messages(x, "mask")		
	}
)


setMethod("project", signature(x="SpatRaster"), 
	function(x, y, method="bilinear", mask=FALSE, filename="", overwrite=FALSE, wopt=list(), ...)  {
		
		method <- ifelse(method == "ngb", "near", method)
		opt <- .runOptions(filename, overwrite, wopt)
		if (inherits(y, "SpatRaster")) {
			#x@ptr <- x@ptr$warp(y@ptr, method, opt)
			x@ptr <- x@ptr$warp(y@ptr, "", method, mask, opt)
		} else {
			if (!is.character(y)) {
				warning("crs should be a character value")
				y <- as.character(crs(y))
			}
			#x@ptr <- x@ptr$warpcrs(y, method, opt)
			x@ptr <- x@ptr$warp(SpatRaster$new(), y, method, mask, opt)
		}
		show_messages(x, "project")
	}
)


setMethod("project", signature(x="SpatVector"), 
	function(x, y, ...)  {
		if (!is.character(y)) {
			y <- crs(y)
		}
		x@ptr <- x@ptr$project(y)
		show_messages(x, "project")
	}
)


setMethod("quantile", signature(x="SpatRaster"), 
	function(x, probs=seq(0, 1, 0.25), na.rm=FALSE, filename="", overwrite=FALSE, wopt=list(), ...) { 
		opt <- .runOptions(filename, overwrite, wopt)
		x@ptr <- x@ptr$quantile(probs, na.rm[1], opt)
		show_messages(x, "quantile")
	}
)


setMethod("rasterize", signature(x="SpatVector", y="SpatRaster"), 
	function(x, y, field=1:nrow(x), background=NA, update=FALSE, touches=is.lines(x), filename="", overwrite=FALSE, wopt=list(), ...) { 
		inverse=FALSE # use "mask" for TRUE
		opt <- .runOptions(filename, overwrite, wopt)
		background <- as.numeric(background[1])
		#if (is.na(background)) background = 0/0 # NAN
		if (is.character(field)) {
			y@ptr <- y@ptr$rasterize(x@ptr, field, 0, background, update[1], touches[1], inverse[1], opt)
		} else if (is.numeric(field)) {
			y@ptr <- y@ptr$rasterize(x@ptr, "", field, background, update[1], touches[1], inverse[1], opt)
		} else {
			stop("field should be character or numeric")
		}
		show_messages(y, "rasterize")
	}
)


setMethod("rotate", signature(x="SpatRaster"), 
	function(x, left=TRUE, filename="", overwrite=FALSE, wopt=list(), ...) { 
		opt <- .runOptions(filename, overwrite, wopt)
		x@ptr <- x@ptr$rotate(left, opt)
		show_messages(x, "rotate")		
	}
)



setMethod("shift", signature(x="SpatRaster"), 
	function(x, dx=0, dy=0, filename="", overwrite=FALSE, wopt=list(), ...) { 
		opt <- .runOptions(filename, overwrite, wopt)
		x@ptr <- x@ptr$shift(dx, dy, opt)
		show_messages(x, "shift")		
	}
)

setMethod("slope", signature(x="SpatRaster"), 
	function(x, neighbors=8, unit="degrees", filename="", overwrite=FALSE, wopt=list(), ...) { 
		opt <- .runOptions(filename, overwrite, wopt)
		stopifnot(neighbors %in% c(4,8))
		stopifnot(unit %in% c("degrees", "radians"))
		x@ptr <- x@ptr$slope(neighbors, unit=="degrees", opt)
		show_messages(x, "slope")		
	}
)




setMethod("t", signature(x="SpatRaster"), 
	function(x) {
		opt <- .runOptions(filename="", overwrite=TRUE, wopt=list())
		x@ptr <- x@ptr$transpose(opt)
		show_messages(x, "t")
	}
)


setMethod("trim", signature(x="SpatRaster"), 
	function(x, padding=0, filename="", overwrite=FALSE, wopt=list(), ...) {
		opt <- .runOptions(filename, overwrite, wopt)
		x@ptr <- x@ptr$trim(padding[1], opt)
		show_messages(x, "rasterize")
	}
)

setMethod("transpose", signature(x="SpatRaster"), 
	function(x, filename="", overwrite=FALSE, wopt=list(), ...) {
		opt <- .runOptions(filename, overwrite, wopt)
		x@ptr <- x@ptr$transpose(opt)
		show_messages(x, "transpose")
	}
)

setMethod("unique", signature(x="SpatRaster", incomparables="ANY"), 
	function(x, incomparables=FALSE, ...) {
		opt <- .getOptions()
		u <- x@ptr$unique(incomparables, opt)
		if (!incomparables) {
			if (!length(u)) return(u)
			u <- do.call(cbind, u)
			colnames(u) = names(x)
		}
		u
	}
)


setMethod("resample", signature(x="SpatRaster", y="SpatRaster"), 
	function(x, y, method="bilinear", filename="", overwrite=FALSE, wopt=list(), ...)  {
		method <- ifelse(method == "ngb", "near", method)
		opt <- .runOptions(filename, overwrite, wopt)
		x@ptr <- x@ptr$warp(y@ptr, "", method, FALSE, opt)
		show_messages(x, "resample")
	}
)

setMethod("summary", signature(object="SpatRaster"), 
	function(object, size=100000, ...)  {
		summary(spatSample(object, size, method="regular", ...))
	}
)


setMethod("stretch", signature(x="SpatRaster"), 
	function(x, minv=0, maxv=255, minq=0, maxq=1, smin=NA, smax=NA, filename="", overwrite=FALSE, wopt=list(), ...) {
		opt <- .runOptions(filename, overwrite, wopt)
		x@ptr <- x@ptr$stretch(minv, maxv, minq, maxq, smin, smax, opt)
		show_messages(x, "stretch")
	}
)



#setMethod("warp", signature(x="SpatRaster", y="SpatRaster"), 
#	function(x, y, method="bilinear", filename="", overwrite=FALSE, wopt=list(), ...)  {
#		opt <- .runOptions(filename, overwrite, wopt)
#		x@ptr <- x@ptr$warp(y@ptr, "", method, opt)
#		show_messages(x, "warp")
#	}
#)

