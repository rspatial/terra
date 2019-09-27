# Author: Robert J. Hijmans
# Date : October 2018
# Version 1.0
# License GPL v3


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
	function(x, filename="", overwrite=FALSE, wopt=list(), ...) {
		opt <- .runOptions(filename, overwrite, wopt)
		x@ptr <- x@ptr$area(opt)
		show_messages(x, "area")
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


setMethod("c", signature(x="SpatRaster"), 
	function(x, ...) {
		dots <- list(...)
		for (i in dots) {
			if (class(i) == "SpatRaster") {
				x@ptr <- x@ptr$combineSources(i@ptr)
			}
		}
		x@ptr$setNames(x@ptr$names)
		show_messages(x, "c")		
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



setMethod("collapse", signature(x="SpatRaster"), 
	function(x, y, filename="", overwrite=FALSE, wopt=list(), ...) { 
		opt <- .runOptions(filename, overwrite, wopt)
		x@ptr <- x@ptr$collapse(y@ptr, opt)
		show_messages(x, "collapse")		
	}
)

setMethod("cover", signature(x="SpatRaster", y="SpatRaster"), 
	function(x, y, value=NA, filename="", overwrite=FALSE, wopt=list(), ...) {
		opt <- .runOptions(filename, overwrite, wopt)
		x@ptr <- x@ptr$cover(y@ptr, value[1], opt)
		show_messages(x, "cover")		
	}
)

setMethod("disaggregate", signature(x="SpatRaster"), 
	function(x, fact, filename="", overwrite=FALSE, wopt=list(), ...) {
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
		v <- x@ptr$freq(bylayer[1])
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
	function(x, mask, inverse=FALSE, updatevalue=NA, filename="", overwrite=FALSE, wopt=list(), ...) { 
		opt <- .runOptions(filename, overwrite,wopt)
		x@ptr <- x@ptr$mask_vector(mask@ptr, inverse[1], NA, updatevalue[1], opt)
		show_messages(x, "mask")		
	}
)



setMethod("project", signature(x="SpatRaster"), 
	function(x, crs, method="bilinear", filename="", overwrite=FALSE, wopt=list(), ...)  {
		opt <- .runOptions(filename, overwrite, wopt)
		if (!is.character(crs)) {
			warning("crs should be a character value")
			crs <- as.character(crs(crs))
		}
		x@ptr <- x@ptr$project(crs, method, opt)
		show_messages(x, "project")
	}
)

setMethod("project", signature(x="SpatVector"), 
	function(x, crs, ...)  {
		if (!is.character(crs)) {
			crs <- crs(x)
		}
		x@ptr <- x@ptr$project(crs)
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
	function(x, y, field=1:nrow(x), background=NA, update=FALSE, filename="", overwrite=FALSE, wopt=list(), ...) { 
		opt <- .runOptions(filename, overwrite, wopt)
		if (is.character(field)) {
			field <- x[[field, drop=TRUE]]
			if (!is.numeric(field)) {
				stop("this is not a numerical variable")
			}
		}
		y@ptr <- y@ptr$rasterize(x@ptr, field, background[1], update[1], opt)
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
		u <- x@ptr$unique(incomparables)
		if (!incomparables) {
			if (!length(u)) return(u)
			u <- do.call(cbind, u)
			colnames(u) = names(x)
		}
		u
	}
)

setMethod("warp", signature(x="SpatRaster", y="SpatRaster"), 
	function(x, y, method="bilinear", filename="", overwrite=FALSE, wopt=list(), ...)  {
		opt <- .runOptions(filename, overwrite, wopt)
		x@ptr <- x@ptr$warp(y@ptr, method, opt)
		show_messages(x, "warp")
	}
)

