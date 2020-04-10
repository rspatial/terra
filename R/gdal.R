

gdal_info <- function(filename, options="", print=TRUE, open_opt="") {
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

