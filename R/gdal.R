
unloadGDALdrivers <- function(x) {
	.removeDriver(x)
}


fileBlocksize <- function(x) {
	v <- x@pntr$getFileBlocksize()
	m <- matrix(v, ncol=2)
	colnames(m) <- c("rows", "cols")
	m
}


clearVSIcache <- function() {
	.clearVSIcache(TRUE)
}

gdalCache <- function(size=NA) {
	vsi <- FALSE # vsi not working
	if (is.null(size) || is.na(size)) {
		.getGDALCacheSizeMB(vsi)
	} else if (size > 0) {
		.setGDALCacheSizeMB(size, vsi)
	}
}


getGDALconfig <- function(option) {
	sapply(option, .gdal_getconfig)
}

setGDALconfig <- function(option, value="") {
	value <- rep_len(value, length(option))
	for (i in 1:length(option)) {
		if (grepl("=", option[i])) {
			opt <- trimws(unlist(strsplit(option[i], "="))[1:2])
			.gdal_setconfig(opt[1], opt[2])
		} else {
			.gdal_setconfig(trimws(option[i]), trimws(value[i]))
		}
	}
}


libVersion <- function(lib="all", parse=FALSE) {
	lib <- tolower(lib)
	if (lib=="gdal") {
		out <- .gdal_version()
	} else if (lib=="proj") {
		out <- proj_version()
	} else if (lib=="geos") {
		out <- .geos_version()
	} else {
		out <- c(gdal=.gdal_version(), proj=proj_version(), geos=.geos_version())
	}
	if (parse) {
		nms <- names(out)
		out <- data.frame(matrix(as.numeric(unlist(strsplit(out, "\\."))), ncol=3, byrow=TRUE), row.names=nms)
		names(out) <- c("major", "minor", "sub")
	}
	out
}




gdal <- function(warn=NA, drivers=FALSE, ...) {
	if (!is.na(warn)) {
		warn <- as.integer(warn)
		stopifnot(warn %in% (1:4))
		.set_gdal_warnings(warn)
	} else if (drivers) {
		x <- .gdaldrivers()
		x <- do.call(cbind, x)
		x <- data.frame(x)
		x[,2] = c(FALSE, TRUE)[as.integer(x[,2])+1]
		x[,3] = c(FALSE, TRUE)[as.integer(x[,3])+1]
		x[,4] = c("read", "read/write", "read/write")[as.integer(x[,4])+1]
		colnames(x) <- c("name", "raster", "vector", "can", "vsi", "long.name")
		x[,5] <- x[,5] == 1
		x <- x[order(x$name), ]
		rownames(x) <- NULL
		x
	} else {
		dots <- list(...)
		if (length(dots) > 0) {
			libVersion(dots[1])
		} else {
			libVersion("gdal")		
		}
	}
}



.describe_sds <- function(x, print=FALSE) {
	x <- .sdinfo(x)
	if (length(x[[1]]) == 1 & length(x[[2]]) == 0) {
		error("gdal (sds)", x[[1]])
	}
	names(x) <- c("name", "var", "desc", "nrow", "ncol", "nlyr")
	m <- do.call(cbind, x)
	m <- data.frame(id=1:nrow(m), m, stringsAsFactors=FALSE)
	ii <- which(colnames(m) %in% c("nrow", "ncol", "nlyr"))
	for (i in ii) m[,i] <- as.integer(m[,i])
	if (print) {
		print(m)
		invisible(m)
	} else {
		m
	}
}



.meta_sds <- function(x, parse=FALSE, ...) {
	if (parse) {
		m <- .parsedsdsmetadata(x)
		m <- do.call(cbind, m)
		if (nrow(m) > 0) {
			m <- data.frame(1:nrow(m), m, stringsAsFactors=FALSE)
		} else {
			m <- data.frame(0[0], m, stringsAsFactors=FALSE)
		}
		for (i in 5:7) m[,i] <- as.integer(m[,i])
		colnames(m) <- c("id", "name", "var", "desc", "nrow", "ncol", "nlyr")
	} else {
		m <- .sdsmetadata(x)
	}
	m
}



setMethod("describe", signature(x="character"),
	function(x, sds=FALSE, meta=FALSE, parse=FALSE, options="", print=FALSE, open_opt="", mdim=FALSE) {

		if (mdim) {
			return(.gdalmdinfo(x, ""[0]))
		}

		#x <- .fullFilename(x[1], FALSE)
		x <- x[1]
		if (meta) {
			if (sds) {
				return(.meta_sds(x, parse))
			} else {
				return(.metadata(x))
			}
		}

		if (sds) {
			return(.describe_sds(x, print=print))
		}
		options <- unique(trimws(options))
		options <- options[options != ""]
		if (length(options) > 0) {
			options <- paste0("-", options)
			options <- gsub("^--", "-", options)
		}
		open_opt <- unique(trimws(open_opt))
		open_opt <- open_opt[open_opt != ""]
		g <- .gdalinfo(x, options, open_opt)

		if (g == "") {
			add <- ifelse(!file.exists(x), "\n", "\nThe file does not exist\n")
			x <- paste0("GDAL cannot open: ", x, add)
		}
		y <- unlist(strsplit(g, "\n"))
		if (print) {
			cat(g, "\n")
			invisible(y)
		} else {
			return(y)
		}
	}
)

setMethod("describe", signature(x="SpatRaster"),
	function(x, source=1, ...) {
		if (!hasValues(x)) return(NULL)
		source <- round(source)
		if ((source < 1) || (source > nsrc(x))) {
			error("describe", "source should be >= 1 and <= nsrc()")
		}
		describe(sources(x)[source])
	}
)
