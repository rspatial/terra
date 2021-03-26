
.terra_environment <- new.env()

 
.create_options <- function() {
	opt <- methods::new("SpatOptions")
	opt@ptr <- SpatOptions$new()
	opt@ptr$tempdir <- normalizePath(tempdir(check = TRUE), winslash="/")
	.terra_environment$options <- opt
}
 
.options_names <- function() {
	c("progress", "tempdir", "memfrac", "datatype", "filetype", "filenames", "overwrite", "todisk", "names", "verbose", "NAflag", "statistics", "steps", "ncopies") 
}

 
.setOptions <- function(x, wopt) {

	nms <- names(wopt)

	g <- which(nms == "gdal")
	if (length(g) > 0) {
		gopt <- unlist(wopt[g])
		wopt <- wopt[-g]
		nms <- nms[-g]
		i <- grep("=", gopt)
		gopt <- gopt[i]
		gopt <- gsub(" ", "", gopt)
		x$gdal_options <- gopt
	}

	s <- nms %in% .options_names()

	if (any(!s)) {
		bad <- paste(nms[!s], collapse=",")
		error("write", "unknown option(s):", bad)
	}

	if (any(s)) {
		nms <- nms[s]
		wopt <- wopt[s]
		i <- which(nms == "names")
		if (length(i) > 0) {
			namevs <- trimws(unlist(strsplit(wopt[[i]], ",")))
			x[["names"]] <- namevs
			wopt <- wopt[-i]
			nms <- nms[-i]
		}

		for (i in seq_along(nms)) {
			x[[nms[i]]] <- wopt[[i]]
		}
		if ("datatype" %in% nms) {
			x$datatype_set = TRUE;
		}
	}

	x
} 
 
 
spatOptions <- function(filename="", overwrite=FALSE, ..., wopt=NULL) {

	w <- list(...)
	wopt <- c(w, wopt)
	
	## work around onLoad problem
	if (is.null(.terra_environment$options)) .create_options()

	ptr <- .terra_environment$options@ptr
	opt <- ptr$deepcopy()

	filename <- .fullFilename(filename, mustExist=FALSE)
	if (!is.null(unlist(wopt))) {
		wopt$filenames <- filename
		wopt$overwrite <- overwrite[1]
		opt <- .setOptions(opt, wopt)
	} else {
		opt$filenames <- filename
		opt$overwrite <- overwrite[1]
	}
	#messages(opt)
	#opt$todisk <- TRUE
	opt
}

#..getOptions <- function() {
#	spatOptions("", TRUE, list())
#}

#..showOptions <- function(opt) {
#	cat("Options for package 'terra'\n")
#	cat("memfrac     :" , opt$memfrac, "\n")
#	cat("tempdir     :" , opt$tempdir, "\n")
#	cat("datatype    :" , opt$def_datatype, "\n")
#	cat("filetype    :" , opt$def_filetype, "\n")
#	cat("progress    :" , opt$progress, "\n")
#	cat("verbose     :" , opt$verbose, "\n")
#	if (opt$todisk) {
#		cat("todisk      :" , opt$todisk, "\n")
#	}
#}

.showOptions <- function(opt) {
	nms <- c("memfrac", "tempdir", "datatype", "progress", "todisk", "verbose") 
	for (n in nms) {
		v <- eval(parse(text=paste0("opt$", n)))
		cat(paste0(substr(paste(n, "         "), 1, 10), ": ", v, "\n"))
	}
}


.default_option_names <- function() {
	c("datatype", "filetype") #, "verbose") 
}


 
terraOptions <- function(...) {
	dots <- list(...)
	opt <- .terra_environment$options@ptr
	if (length(dots) == 0) {
		.showOptions(opt)
	} else {
		nms <- names(dots)
		d <- nms %in% .default_option_names()
		dnms <- paste0("def_", nms)
		for (i in 1:length(nms)) {
			if (d[i]) {
				opt[[ dnms[i] ]] <- dots[[ i ]]
			} else {
				opt[[ nms[i] ]] <- dots[[ i ]]
			}
		}
		if ("memfrac" %in% nms) {
			if (dots$memfrac > 0.9) {
				warn("terraOptions", "memfrac > 0.9")
			}
		}
		.terra_environment$options@ptr <- opt
	}
}

