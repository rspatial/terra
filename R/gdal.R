

gdal_info <- function(filename, options="", print=TRUE, open_opt="", ...) {
	if (is.null(.terra_environment$options)) .init()

	options <- unique(trimws(options))
	options <- options[options != ""]
	if (length(options) > 0) {
		options <- paste0("-", options)
		options <- gsub("^--", "-", options)
	}
	open_opt <- unique(trimws(open_opt))
	open_opt <- open_opt[open_opt != ""]
	x <- .gdalinfo(filename, options, open_opt)
	if (print) {
		if (x == "") {
			add <- ifelse(file.exists(filename), "\n", "\nThe file does not exist\n")
			x <- paste0("GDAL cannot open: ", filename, add)
		}
		cat(x, "\n")
		invisible(x)
	} else {
		return(x)
	}
}


sds_info <- function(filename, ...) {
	if (is.null(.terra_environment$options)) .init()
	x <- .sdinfo(filename)
	if (length(x[[1]]) == 1 & length(x[[2]]) == 0) {
		stop(x[[1]])
	}
	m <- do.call(cbind, x)
	m <- data.frame(1:nrow(m), m, stringsAsFactors=FALSE)
	colnames(m) <- c("id", "name", "desc", "nrow", "ncol", "nlyr")
	for (i in 4:6) m[,i] <- as.integer(m[,i])
	m
}

