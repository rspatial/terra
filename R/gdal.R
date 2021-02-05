

gdal_version <- function() {
	gdal()
}


gdal <- function(warn=NA, drivers=FALSE) {
	if (!is.na(warn)) {
		warn <- as.integer(warn)
		stopifnot(warn %in% c(1:4))
		.set_gdal_warnings(warn)
	} else if (drivers) {
		x <- .gdaldrivers()
		x <- do.call(cbind, x)
		x[,2] = c("vector", "raster")[as.integer(x[,2])+1]
		x[,3] = c("read", "read/copy-write", "read/write")[as.integer(x[,3])+1]
		colnames(x) <- c("name", "type", "can", "long.name")
		x <- data.frame(x)
		x <- x[order(x$type, x$name), ]
		rownames(x) <- NULL
		x
	} else {
		.gdalversion()
	}
}



.describe_sds <- function(x, print=FALSE) {
	x <- .sdinfo(x)
	if (length(x[[1]]) == 1 & length(x[[2]]) == 0) {
		error("gdal (sds)", "not working for: ", x[[1]])
	}
	m <- do.call(cbind, x)
	m <- data.frame(1:nrow(m), m, stringsAsFactors=FALSE)
	colnames(m) <- c("id", "name", "var", "desc", "nrow", "ncol", "nlyr")
	for (i in 5:7) m[,i] <- as.integer(m[,i])
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

