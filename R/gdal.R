

gdal <- function(warn = NA) {
	if (!is.na(warn)) {
		warn <- as.integer(warn)
		stopifnot(warn %in% c(1:4))
		.set_gdal_warnings(warn)
	} else {
		.gdalversion()
	}
}



.describe_sds <- function(x, print=FALSE, ...) {
	x <- .sdinfo(x)
	if (length(x[[1]]) == 1 & length(x[[2]]) == 0) {
		stop(x[[1]])
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

setMethod("desc", signature(x="character"), 
	function(x, ...) {
		stop("depracated function. Use 'describe'")
	}
)


setMethod("describe", signature(x="character"), 
	function(x, sds=FALSE, meta=FALSE, parse=FALSE, options="", print=FALSE, open_opt="", ...) {

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
		x <- .gdalinfo(x, options, open_opt)

		if (x == "") {
			add <- ifelse(file.exists(filename), "\n", "\nThe file does not exist\n")
			x <- paste0("GDAL cannot open: ", x, add)
		}
		y <- unlist(strsplit(x, "\n"))
		if (print) {
			cat(x, "\n")
			invisible(y)
		} else {
			return(y)
		}
	}
)

