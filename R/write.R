
setMethod("writeStart", signature(x="SpatRaster", filename="character"), 
	function(x, filename="", overwrite=FALSE, n=4, ...) {
		opt <- spatOptions(filename, overwrite, ncopies=n, ...)
		ok <- x@ptr$writeStart(opt)
		messages(x, "writeStart")
		b <- x@ptr$getBlockSize(n, opt$memfrac)
		b$row <- b$row + 1
		b
	}
)


setMethod("writeStop", signature(x="SpatRaster"), 
	function(x) {
		success <- x@ptr$writeStop()
		messages(x, "writeStop")
		f <- sources(x)$source
		if (f != "") {
			x <- rast(f)
		}
		return(x)
	} 
)

setMethod("writeValues", signature(x="SpatRaster", v="vector"), 
	function(x, v, start, nrows) {
		#wstart <- start[1]-1
		#nrows <- start[2]
		#if (is.na(nrows)) {
		#	nrows <- length(v) / (ncol(x) * nlyr(x))
		#}
		success <- x@ptr$writeValues(v, start-1, nrows, 0, ncol(x))
		messages(x, "writeValues")
		invisible(success)
	}
)


setMethod("writeRaster", signature(x="SpatRaster", filename="character"), 
function(x, filename="", overwrite=FALSE, ...) {
	filename <- trimws(filename)
	stopifnot(filename != "")
	if (tools::file_ext(filename) %in% c("nc", "cdf") || isTRUE(list(...)$filetype=="netCDF")) {
		warn("consider writeCDF to write ncdf files")
	}
	opt <- spatOptions(filename, overwrite, ...)
	x@ptr <- x@ptr$writeRaster(opt)
	x <- messages(x, "writeRaster")
	invisible(rast(filename))
}
)


setMethod("writeVector", signature(x="SpatVector", filename="character"), 
function(x, filename, filetype="ESRI Shapefile", overwrite=FALSE, options=NULL) {
	filename <- trimws(filename)
	if (filename == "") {
		error("writeVector", "provide a filename")
	}
	
	lyrname <- tools::file_path_sans_ext(basename(filename))
	if (is.null(options)) { options <- ""[0] }
	success <- x@ptr$write(filename, lyrname, filetype, overwrite[1], options)
	messages(x, "writeVector")
	invisible(TRUE)
}
)

