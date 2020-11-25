

gdal_version <- function() {
	#if (is.null(.terra_environment$options)) .init()
	.gdalversion()
}

gdal_warnings <- function(level = 2) {
	level <- as.integer(level)
	stopifnot(level %in% c(1:4))
	.set_gdal_warnings(level - 1)
}


describe <- function(x, options="", print=TRUE, open_opt="", ...) {
	#if (is.null(.terra_environment$options)) .init()

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


describe_sds <- function(x, ...) {
	#if (is.null(.terra_environment$options)) .init()
	x <- .sdinfo(x)
	if (length(x[[1]]) == 1 & length(x[[2]]) == 0) {
		stop(x[[1]])
	}
	m <- do.call(cbind, x)
	m <- data.frame(1:nrow(m), m, stringsAsFactors=FALSE)
	colnames(m) <- c("id", "name", "var", "desc", "nrow", "ncol", "nlyr")
	for (i in 5:7) m[,i] <- as.integer(m[,i])
	m
}



setMethod("meta", signature(x="character"), 
	function(x, ...) {
		.metadata(x)
	}
)

meta_sds <- function(x, parse=FALSE, ...) {
	if (parse) {
		m <- .parsedsdsmetadata(x)
		m <- do.call(cbind, m)
		m <- data.frame(1:nrow(m), m, stringsAsFactors=FALSE)
		colnames(m) <- c("id", "name", "var", "desc", "nrow", "ncol", "nlyr")
		for (i in 5:7) m[,i] <- as.integer(m[,i])	
	} else {
		m <- .sdsmetadata(x)
	}
	m
}

