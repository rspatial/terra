# Author: Robert J. Hijmans
# Date : October 2018
# Version 1.0
# License GPL v3


setMethod("c", signature(x="SpatRaster"), 
	function(x, ...) {
		for (i in list(...)) {
			if (class(i) == "SpatRaster") {
				x@ptr <- x@ptr$combineSources(i@ptr)
			}
		}
		x@ptr$setNames(x@ptr$names)
		show_messages(x, "c")		
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


setMethod("area", signature(x="SpatRaster"), 
	function(x, filename="", overwrite=FALSE, wopt=list(), ...) {
		opt <- .runOptions(filename, overwrite, wopt)
		x@ptr <- x@ptr$area(opt)
		show_messages(x, "area")
	}
)


setMethod("clamp", signature(x="SpatRaster"), 
	function(x, lower=-Inf, upper=Inf, values=TRUE, filename="", overwrite=FALSE, wopt=list(), ...) {
		opt <- .runOptions(filename, overwrite,wopt)
		x@ptr <- x@ptr$clamp(lower, upper, values[1], opt)
		show_messages(x, "clamp")
	}
)

setMethod("crop", signature(x="SpatRaster", y="ANY"), 
	function(x, y, snap="near", filename="", overwrite=FALSE, wopt=list(), ...) {
		if (!inherits(y, "SpatExtent")) {
			y <- try(ext(y))
			if (class(y) == "try-error") { stop("cannot get an extent from y") }
		}
		opt <- .runOptions(filename, overwrite, wopt)
		x@ptr <- x@ptr$crop(y@ptr, snap[1], opt)
		show_messages(x, "crop")		
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

setMethod("gridDistance", signature(x="SpatRaster"), 
	function(x, filename="", overwrite=FALSE, wopt=list(), ...) {
		opt <- .runOptions(filename, overwrite,wopt)
		x@ptr <- x@ptr$gridDistance(opt)
		show_messages(x, "gridDistance")
	}
)

setMethod("flip", signature(x="SpatRaster"), 
	function(x, vertical=TRUE, filename="", overwrite=FALSE, wopt=list(), ...) {
		opt <- .runOptions(filename, overwrite, wopt)
		x@ptr <- x@ptr$flip(vertical, opt)
		show_messages(x, "flip")
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

setMethod("merge", signature(x="SpatRaster", y="SpatRaster"), 
	function(x, y, filename="", overwrite=FALSE, wopt=list(), ...) { 
		opt <- .runOptions(filename, overwrite, wopt)
		x@ptr <- x@ptr$merge(y, opt)
		show_messages(x, "merge")		
	}
)


setMethod("project", signature(x="SpatRaster"), 
	function(x, crs, method="bilinear", filename="", overwrite=FALSE, wopt=list(), ...)  {
		opt <- .runOptions(filename, overwrite, wopt)
		x@ptr <- x@ptr$project(crs, method, opt)
		show_messages(x, "project")
	}
)

setMethod("project", signature(x="SpatVector"), 
	function(x, crs, ...)  {
		x@ptr <- x@ptr$project(crs)
		show_messages(x, "project")
	}
)

setMethod("rasterize", signature(x="SpatVector", y="SpatRaster"), 
	function(x, y, background=NA, filename="", overwrite=FALSE, wopt=list(), ...) { 
		opt <- .runOptions(filename, overwrite, wopt)
		y@ptr <- y@ptr$rasterize(x@ptr, background[1], opt)
		show_messages(y, "rasterize")
	}
)


setMethod("reclassify", signature(x="SpatRaster", rcl="ANY"), 
function(x, rcl, include.lowest=FALSE, right=TRUE, othersNA=FALSE, filename="", overwrite=FALSE, wopt=list(), ...) {
	
	if ( is.null(dim(rcl)) ) { 
		stopifnot((length(rcl) %% 3 == 0))
		rcl <- matrix(rcl, ncol=3, byrow=TRUE) 
	} else if (is.data.frame(rcl)) {
		rcl <- as.matrix(rcl)
	}

	right <- ifelse(is.na(right), 2, ifelse(right, 0, 1))
	include.lowest <- as.logical(include.lowest[1])

	opt <- .runOptions(filename, overwrite, wopt)
    x@ptr <- x@ptr$rcppReclassify(rcl, right, include.lowest, othersNA, opt)
	show_messages(x, "reclassify")	
}
)



setMethod("rotate", signature(x="SpatRaster"), 
	function(x, left=TRUE, filename="", overwrite=FALSE, wopt=list(), ...) { 
		opt <- .runOptions(filename, overwrite, wopt)
		x@ptr <- x@ptr$rotate(left, opt)
		show_messages(x, "rotate")		
	}
)


setMethod("sampleRegular", signature(x="SpatRaster", size="numeric"), 
	function(x, size, ...) { 
		size <- max(1, min(size(x), size))
		x@ptr <- x@ptr$sampleRegular(size)
		show_messages(x, "sampleRegular")		
	}
)


setMethod("shift", signature(x="SpatRaster"), 
	function(x, dx=0, dy=0, filename="", overwrite=FALSE, wopt=list(), ...) { 
		opt <- .runOptions(filename, overwrite, wopt)
		x@ptr <- x@ptr$shift(dx, dy, opt)
		show_messages(x, "shift")		
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

setMethod("t", signature(x="SpatRaster"), 
	function(x) {
		opt <- .runOptions(filename="", overwrite=TRUE, wopt=list())
		x@ptr <- x@ptr$transpose(opt)
		show_messages(x, "t")
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
