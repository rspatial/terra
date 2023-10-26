

fileBlocksize <- function(x) {
	v <- x@cpp$getFileBlocksize()
	m <- matrix(v, ncol=2)
	colnames(m) <- c("rows", "cols")
	m
}


gdalCache <- function(size=NA) {
	if (is.na(size)) {
		.getGDALCacheSizeMB()
	} else {
		if (size > 0) {
			.setGDALCacheSizeMB(size)
		}
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




gdal <- function(warn=NA, drivers=FALSE, lib="gdal") {
	if (!is.na(warn)) {
		warn <- as.integer(warn)
		stopifnot(warn %in% c(1:4))
		.set_gdal_warnings(warn)
	} else if (drivers) {
		x <- .gdaldrivers()
		x <- do.call(cbind, x)
		x[,2] = c("vector", "raster")[as.integer(x[,2])+1]
		#x[,3] = c("read", "read/copy-write", "read/write")[as.integer(x[,3])+1]
		x[,3] = c("read", "read/write", "read/write")[as.integer(x[,3])+1]
		colnames(x) <- c("name", "type", "can", "vsi", "long.name")
		x <- data.frame(x)
		x[,4] <- x[,4] == 1
		x <- x[order(x$type, x$name), ]
		rownames(x) <- NULL
		x
	} else {
		lib <- tolower(lib)
		if (lib=="gdal") {
			.gdal_version()
		} else if (lib=="proj") {
			proj_version()
		} else if (lib=="geos") {
			.geos_version()
		} else {
			c(gdal=.gdal_version(), proj=proj_version(), geos=.geos_version())
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
	function(x, sds=FALSE, meta=FALSE, parse=FALSE, options="", print=FALSE, open_opt="") {

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

