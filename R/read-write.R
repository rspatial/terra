# Author: Robert J. Hijmans
# Date : December 2017
# Version 1.0
# License GPL v3



setMethod("readStart", signature(x="SpatRaster"), 
	function(x, ...) {
		success <- x@ptr$readStart()
		show_messages(x, "readStart")		
		invisible(success)
	}
)

setMethod("readStart", signature(x="SpatStack"), 
	function(x, ...) {
		success <- x@ptr$readStart()
		show_messages(x, "readStart")		
		invisible(success)
	}
)

#setMethod("readStart", signature(x="SpatStack"), 
#	function(x, ...) {
#		nsd <- nsds(x)
#		for (i in 1:nsd) {
#			y <- x[i]
#			success <- readStart(y)
#			x[i] <- y
#		}
#		show_messages(x, "readStart")		
#		invisible(success)
#	}
#)


setMethod("readStop", signature(x="SpatRaster"), 
	function(x) {
		success <- x@ptr$readStop()
		show_messages(x, "readStop")		
		invisible(success)
	}
)

setMethod("readStop", signature(x="SpatStack"), 
	function(x) {
		success <- x@ptr$readStop()
		show_messages(x, "readStop")		
		invisible(success)
	}
)


setMethod("writeStart", signature(x="SpatRaster", filename="character"), 
	function(x, filename="", overwrite=FALSE, wopt=list(), ...) {
		opt <- .runOptions(filename, overwrite, wopt)
		ok <- x@ptr$writeStart(opt)
		show_messages(x, "writeStart")		
		b <- x@ptr$getBlockSize(4, opt$memfrac)
		b$row <- b$row + 1
		b		
	}
)


setMethod("writeStop", signature(x="SpatRaster"), 
	function(x) {
		success <- x@ptr$writeStop()
		show_messages(x, "writeStop")
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
		show_messages(x, "writeValues")
		invisible(success)
	}
)


setMethod("writeRaster", signature(x="SpatRaster", filename="character"), 
function(x, filename="", overwrite=FALSE, wopt=list(), ...) {
	opt <- .runOptions(filename, overwrite, wopt)
	x@ptr <- x@ptr$writeRaster(opt)
	show_messages(x, "writeRaster")
}
)


setMethod("writeVector", signature(x="SpatVector", filename="character"), 
function(x, filename, overwrite=FALSE, ...) {
	filename <- trimws(filename)
	if (filename == "") {
		stop("provide a filename")
	}
	lyrname <- gsub(".shp", "", basename(filename))
	success <- x@ptr$write(filename, lyrname, "ESRI Shapefile", overwrite[1])
	show_messages(x, "writeVector")
	invisible(TRUE)
}
)

