# Author: Robert J. Hijmans
# Date : December 2017
# Version 1.0
# License GPL v3

readAll <- function(x) {
	ok <- x@ptr$readAll()
	x <- messages(x)
	invisible(ok)
}

setMethod("readStart", signature(x="SpatRaster"), 
	function(x, ...) {
		success <- x@ptr$readStart()
		messages(x, "readStart")		
		if (!success) error("readStart,SpatRaster", "cannot open file for reading")
		invisible(success)
	}
)

setMethod("readStart", signature(x="SpatRasterDataset"), 
	function(x, ...) {
		success <- x@ptr$readStart()
		messages(x, "readStart")		
		if (!success) error("readStart,SpatRasterDataset", "cannot open file for reading")
		invisible(success)
	}
)

#setMethod("readStart", signature(x="SpatRasterDataset"), 
#	function(x, ...) {
#		nsd <- length(x)
#		for (i in 1:nsd) {
#			y <- x[i]
#			success <- readStart(y)
#			x[i] <- y
#		}
#		messages(x, "readStart")		
#		invisible(success)
#	}
#)


setMethod("readStop", signature(x="SpatRaster"), 
	function(x) {
		success <- x@ptr$readStop()
		messages(x, "readStop")		
		invisible(success)
	}
)

setMethod("readStop", signature(x="SpatRasterDataset"), 
	function(x) {
		success <- x@ptr$readStop()
		messages(x, "readStop")		
		invisible(success)
	}
)


setMethod("writeStart", signature(x="SpatRaster", filename="character"), 
	function(x, filename="", overwrite=FALSE, wopt=list(), ...) {
		opt <- spatOptions(filename, overwrite, wopt)
		ok <- x@ptr$writeStart(opt)
		messages(x, "writeStart")		
		b <- x@ptr$getBlockSize(4, opt$memfrac)
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
	function(x, v, start, nrows, ...) {
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
function(x, filename="", overwrite=FALSE, wopt=list(), ...) {
	filename <- trimws(filename)
	stopifnot(filename != "")
	if (tools::file_ext(filename) %in% c("nc", "cdf") || isTRUE(wopt$filetype=="netCDF")) {
		warn("writeRaster", "call writeCDF directly")
		return ( writeCDF(x, filename=filename, overwrite=overwrite, ...) )
	}
	opt <- spatOptions(filename, overwrite, wopt)
	x@ptr <- x@ptr$writeRaster(opt)
	x <- messages(x, "writeRaster")
	invisible(rast(filename))
}
)


setMethod("writeVector", signature(x="SpatVector", filename="character"), 
function(x, filename, overwrite=FALSE, ...) {
	filename <- trimws(filename)
	if (filename == "") {
		error("writeVector", "provide a filename")
	}
	lyrname <- gsub(".shp", "", basename(filename))
	success <- x@ptr$write(filename, lyrname, "ESRI Shapefile", overwrite[1])
	messages(x, "writeVector")
	invisible(TRUE)
}
)

