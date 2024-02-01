# Author: Robert J. Hijmans
# Date :  October 2017
# Version 0.9
# License GPL v3

setMethod("ext", signature(x="SpatExtent"),
	function(x){
		x@ptr <- x@ptr$deepcopy()
		x
	}
)

setMethod("ext", signature(x="SpatRasterCollection"),
	function(x){
		e <- methods::new("SpatExtent")
		e@ptr <- x@ptr$extent()
		e
	}
)


setMethod("ext", signature(x="sf"),
	function(x){
		sfi <- attr(x, "sf_column")
		geom <- x[[sfi]]
		e <- attr(geom, "bbox")
		ext(e[c(1,3,2,4)])
	}
)

setMethod("ext", signature(x="missing"),
	function(x){
		e <- methods::new("SpatExtent")
		e@ptr <- SpatExtent$new()
		e
	}
)

setMethod("ext", signature(x="numeric"),
	function(x, ..., xy=FALSE){
		dots <- as.vector(unlist(list(...)))
		x <- c(x, dots)
		n <- length(x)
		if (n != 4) {
			error("ext", "expected four numbers")
		}
		names(x) <- NULL
		e <- methods::new("SpatExtent")
		if (xy) {
			e@ptr <- SpatExtent$new(x[1], x[3], x[2], x[4])
		} else {
			e@ptr <- SpatExtent$new(x[1], x[2], x[3], x[4])
		}
		if (!e@ptr$valid) {
			error("ext", "invalid extent")
		}
		e
	}
)

setMethod("ext", signature(x="matrix"),
	function(x){
		if ((ncol(x) == 2) && (nrow(x) > 2)) {
			x <- apply(x, 2, range, na.rm=TRUE)
		}
		ext(as.vector(x))
	}
)

setMethod("ext", signature(x="bbox"),
	function(x){
		ext(x[c(1,3,2,4)])
	}
)


setMethod("ext", signature(x="SpatRaster"),
	function(x, cells=NULL){
		if (!is.null(cells)) {
			cells <- stats::na.omit(unique(round(cells)))
			cells <- cells[cells > 0 & cells <= ncell(x)]
			if (length(cells) < 1) {
				stop("no valid cells")
			}
			r <- res(x)
			dx <- r[1] * c(-0.5, 0.5)
			dy <- r[2] * c(-0.5, 0.5)
			ext(range(xFromCell(x, cells)) + dx, range(yFromCell(x, cells)) + dy)
		} else {
			e <- methods::new("SpatExtent")
			e@ptr <- x@ptr$extent
			return(e)
		}
	}
)



setMethod("ext", signature(x="SpatRasterDataset"),
	function(x){
		e <- methods::new("SpatExtent")
		e@ptr <- x@ptr$ext()
		return(e)
	}
)



setMethod("ext<-", signature("SpatRaster", "SpatExtent"),
	function(x, value) {
		x@ptr$extent <- value@ptr
		messages(x, "ext<-")
	}
)


setMethod("ext<-", signature("SpatRaster", "numeric"),
	function(x, value) {
		e <- ext(value)
		x@ptr <- x@ptr$deepcopy()
		x@ptr$extent <- e@ptr
		messages(x, "ext<-")
	}
)

setMethod("set.ext", signature("SpatRaster"),
	function(x, value) {
		e <- ext(value)
		x@ptr$extent <- e@ptr
		messages(x, "set_ext")
		invisible(TRUE)
	}
)


setMethod("ext", signature(x="SpatVector"),
	function(x) {
		e <- methods::new("SpatExtent")
		e@ptr <- x@ptr$extent()
		e
	}
)

setMethod("ext", signature(x="SpatVectorProxy"),
	function(x) {
		e <- methods::new("SpatExtent")
		e@ptr <- x@ptr$v$extent()
		e
	}
)


setMethod("ext", signature(x="Extent"),
	function(x) {
		ext(as.vector(x))
	}
)

setMethod("ext", signature(x="Raster"),
	function(x) {
		ext(x@extent)
	}
)

setMethod("ext", signature(x="Spatial"),
	function(x) {
		ext(as.vector(t(x@bbox)))
	}
)



setMethod("xmin", signature(x="SpatExtent"),
	function(x){
		x@ptr$vector[1]
	}
)
setMethod("xmax", signature(x="SpatExtent"),
	function(x){
		x@ptr$vector[2]
	}
)
setMethod("ymin", signature(x="SpatExtent"),
	function(x){
		x@ptr$vector[3]
	}
)
setMethod("ymax", signature(x="SpatExtent"),
	function(x){
		x@ptr$vector[4]
	}
)


setMethod("xmin<-", signature("SpatExtent", "numeric"),
	function(x, value){
		v <- as.vector(x)
		v[1] <- value
		ext(v)
	}
)
setMethod("xmax<-", signature("SpatExtent", "numeric"),
	function(x, value){
		v <- as.vector(x)
		v[2] <- value
		ext(v)
	}
)
setMethod("ymin<-", signature("SpatExtent", "numeric"),
	function(x, value){
		v <- as.vector(x)
		v[3] <- value
		ext(v)
	}
)
setMethod("ymax<-", signature("SpatExtent", "numeric"),
	function(x, value){
		v <- as.vector(x)
		v[4] <- value
		ext(v)
	}
)


setMethod("xmin", signature(x="SpatRaster"),
	function(x){
		xmin(ext(x))
	}
)
setMethod("xmax", signature(x="SpatRaster"),
	function(x){
		xmax(ext(x))
	}
)
setMethod("ymin", signature(x="SpatRaster"),
	function(x){
		ymin(ext(x))
	}
)
setMethod("ymax", signature(x="SpatRaster"),
	function(x){
		ymax(ext(x))
	}
)


setMethod("xmin<-", signature("SpatRaster", "numeric"),
	function(x, value){
		v <- as.vector(ext(x))
		v[1] <- value
		x@ptr <- x@ptr$deepcopy()
		ext(x) <- ext(v)
		x
	}
)


setMethod("xmax<-", signature("SpatRaster", "numeric"),
	function(x, value){
		v <- as.vector(ext(x))
		v[2] <- value
		x@ptr <- x@ptr$deepcopy()
		ext(x) <- ext(v)
		x
	}
)
setMethod("ymin<-", signature("SpatRaster", "numeric"),
	function(x, value){
		v <- as.vector(ext(x))
		v[3] <- value
		x@ptr <- x@ptr$deepcopy()
		ext(x) <- ext(v)
		x
	}
)

setMethod("ymax<-", signature("SpatRaster", "numeric"),
	function(x, value){
		v <- as.vector(ext(x))
		v[4] <- value
		x@ptr <- x@ptr$deepcopy()
		ext(x) <- ext(v)
		x
	}
)



setMethod("xmin", signature(x="SpatVector"),
	function(x){
		xmin(ext(x))
	}
)
setMethod("xmax", signature(x="SpatVector"),
	function(x){
		xmax(ext(x))
	}
)
setMethod("ymin", signature(x="SpatVector"),
	function(x){
		ymin(ext(x))
	}
)
setMethod("ymax", signature(x="SpatVector"),
	function(x){
		ymax(ext(x))
	}
)



setMethod("$", "SpatExtent",
	function(x, name) {
		as.vector(x)[name]
	}
)


setMethod("$<-", "SpatExtent",
	function(x, name, value) {
		e <- as.vector(x)
		e[name] <- value
		ext(e)
	}
)



setMethod("[", c("SpatExtent", "missing", "missing"),
	function(x, i, j) {
		as.vector(x)
	}
)

setMethod("[", c("SpatExtent", "numeric", "missing"),
	function(x, i, j) {
		x <- as.vector(x)
		x[i]
	}
)

setReplaceMethod("[", c("SpatExtent", "numeric", "missing"),
	function(x, i, j, value) {
		e <- as.vector(x)
		stopifnot(all(i %in% 1:4))
		e[i] <- value
		ext(e)
	}
)


setMethod("is.valid", signature(x="SpatExtent"),
	function(x) {
		x@ptr$valid
		#x@ptr$valid_notempty
	}
)


setMethod("is.empty", signature(x="SpatExtent"),
	function(x) {
		x@ptr$empty
	}
)
