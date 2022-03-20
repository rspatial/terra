

# not exported
if (!isGeneric("blockSize")) {setGeneric("blockSize", function(x, ...) standardGeneric("blockSize"))}
setMethod("blockSize", signature(x="SpatRaster"), 
	function(x, n) {
		opt <- spatOptions("", FALSE, ncopies=n)
		b <- x@ptr$getBlockSizeR(n, opt$memfrac)
		b$row <- b$row + 1
		b
	}
)


setMethod("writeStart", signature(x="SpatRaster", filename="character"), 
	function(x, filename="", overwrite=FALSE, n=4, ...) {
		filename <- enc2utf8(filename)
		opt <- spatOptions(filename, overwrite, ncopies=n, ...)
		ok <- x@ptr$writeStart(opt)
		messages(x, "writeStart")
		b <- x@ptr$getBlockSizeWrite()
		b$row <- b$row + 1
		b
	}
)


setMethod("writeStop", signature(x="SpatRaster"), 
	function(x) {
		success <- x@ptr$writeStop()
		messages(x, "writeStop")
		f <- sources(x)
		if (f != "") {
			x <- rast(f)
		}
		return(x)
	} 
)

setMethod("writeValues", signature(x="SpatRaster", v="vector"), 
	function(x, v, start, nrows) {
		success <- x@ptr$writeValues(v, start-1, nrows)
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
	filename <- enc2utf8(filename)
	opt <- spatOptions(filename, overwrite, ...)
	x@ptr <- x@ptr$writeRaster(opt)
	x <- messages(x, "writeRaster")
	invisible(rast(filename))
}
)


get_filetype <- function(filename) {
	ext <- tolower(tools::file_ext(filename))
	if (ext == "shp" || ext == "") {
		"ESRI Shapefile"
	} else if (ext == "gpkg") {
		"GPKG"
	} else if (ext == "gml") {
		"GML"
	} else if (ext == "json") {
		"GeoJSON"
	} else if (ext == "cdf") {
		"netCDF"
	} else if (ext == "svg") {
		"SVG"
	} else if (ext == "kml") {
		"KML"
	} else if (ext == "vct") {
		"Idrisi"
	} else {
		error("writeVector", "cannot guess filetype from filename")
	}
}

setMethod("writeVector", signature(x="SpatVector", filename="character"), 
function(x, filename, filetype=NULL, layer=NULL, insert=FALSE, overwrite=FALSE, options="ENCODING=UTF-8") {
	filename <- trimws(filename)
	filename <- enc2utf8(filename)
	if (filename == "") {
		error("writeVector", "provide a filename")
	}
	if (is.null(filetype)) {
		filetype <- get_filetype(filename)
	}
	if (is.null(layer)) layer <- tools::file_path_sans_ext(basename(filename))
	layer <- trimws(layer)
	
	if (is.null(options)) { options <- ""[0] }

	if (filetype == "ESRI Shapefile") {
		nms <- names(x)
		i <- nchar(nms) > 10
		if (any(i)) {
			nms[i] <- substr(nms[i], 1, 10)
			testnms <- make.unique(nms, sep="")
			if (!all(testnms == nms)) {
				i <- which(i)
				newnms <- substr(nms[i], 1, 9)
				newnms <- make.unique(newnms, sep="")
				j <- which(nchar(newnms) == 9)
				newnms[j] <- paste0(newnms[j], "0")
				nms[i] <- newnms
			}
			x@ptr <- x@ptr$deepcopy()
			names(x) <- nms
		}
	}
	success <- x@ptr$write(filename, layer, filetype, insert[1], overwrite[1], options)
	messages(x, "writeVector")
	invisible(TRUE)
}
)



# setMethod("writeVector", signature(x="SpatVectorProxy", filename="character"), 
# function(x, filename, filetype=NULL, layer=NULL, insert=FALSE, overwrite=FALSE, options="ENCODING=UTF-8") {
	# filename <- trimws(filename)
	# filename <- enc2utf8(filename)
	# if (filename == "") {
		# error("writeVector", "provide a filename")
	# }
	# if (is.null(filetype)) {
		# filetype <- get_filetype(filename)
	# }
	# if (is.null(layer)) layer <- tools::file_path_sans_ext(basename(filename))
	# layer <- trimws(layer)
	
	# if (is.null(options)) { options <- ""[0] }

	# if (filetype == "ESRI Shapefile") {
		# nms <- names(x)
		# i <- nchar(nms) > 10
		# if (any(i)) {
			# nms[i] <- substr(nms[i], 1, 10)
			# testnms <- make.unique(nms, sep="")
			# if (!all(testnms == nms)) {
				# i <- which(i)
				# newnms <- substr(nms[i], 1, 9)
				# newnms <- make.unique(newnms, sep="")
				# j <- which(nchar(newnms) == 9)
				# newnms[j] <- paste0(newnms[j], "0")
				# nms[i] <- newnms
			# }
			# x@ptr <- x@ptr$deepcopy()
			# names(x) <- nms
		# }
	# }
	# success <- x@ptr$write_proxy(filename, layer, filetype, insert[1], FALSE, overwrite[1], options)
	# messages(x, "writeVector")
	# invisible(TRUE)
# }
# )

