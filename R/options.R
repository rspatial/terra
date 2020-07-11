
.terra_environment <- new.env()

 
.create_options <- function() {
	opt <- methods::new("SpatOptions")
	opt@ptr <- SpatOptions$new()
	opt@ptr$tempdir <- normalizePath(tempdir(check = TRUE), winslash="/")
	.terra_environment$options <- opt
}
 
.options_names <- function() {
	c("progress", "tempdir", "memfrac", "datatype", "filetype", "filename", "overwrite", "todisk", "names", "verbose", "NAflag") 
}

 
.setOptions <- function(x, opt) {
	nms <- names(opt)
	
	g <- which(nms == "gdal")
	if (length(g) > 0) {
		gopt <- unlist(opt[g])
		opt <- opt[-g]
		nms <- nms[-g]
		i <- grep("=", gopt)
		gopt <- gopt[i]
		gopt <- gsub(" ", "", gopt)
		x$gdal_options <- gopt
	}
	
	s <- nms %in% .options_names()
	
	if (any(!s)) {
		bad <- paste(nms[!s], collapse=",")
		warning(paste("unknown write option(s):", bad), call. = FALSE)
	}
		
	if (any(s)) {
		nms <- nms[s]
		opt <- opt[s]
		i <- which(nms == "names")	
		if (length(i) > 0) {
			namevs <- trimws(unlist(strsplit(opt[[i]], ",")))
			x[["names"]] <- namevs
			opt <- opt[-i]
			nms <- nms[-i]
		}
		
		for (i in seq_along(nms)) {
			x[[nms[i]]] <- opt[[i]]
		}
	}
	x
} 
 
.runOptions <- function(filename="", overwrite=FALSE, wopt=list()) {
	if (!is.list(wopt)) {
		stop("wopt must be a list")
	}
	
	## work around onLoad problem
	if (is.null(.terra_environment$options)) .create_options()
	
	ptr <- .terra_environment$options@ptr
	opt <- ptr$deepcopy(ptr)
	
	filename <- .fullFilename(filename[1], mustExist=FALSE)
	if (!is.null(unlist(wopt))) {
		wopt$filename <- filename
		wopt$overwrite <- overwrite[1]
		opt <- .setOptions(opt, wopt)
	} else {
		opt$filename <- filename
		opt$overwrite <- overwrite[1]
	}
	#show_messages(opt)
	#opt$todisk <- TRUE
	opt
}

.getOptions <- function() {
	.runOptions("", TRUE, list())
}

..showOptions <- function(opt) {
	cat("Options for package 'terra'\n")
	cat("memfrac     :" , opt$memfrac, "\n")
	cat("tempdir     :" , opt$tempdir, "\n")
	cat("datatype    :" , opt$def_datatype, "\n")
	cat("filetype    :" , opt$def_filetype, "\n")
	cat("progress    :" , opt$progress, "\n")
	cat("verbose     :" , opt$verbose, "\n")
	if (opt$todisk) {
		cat("todisk      :" , opt$todisk, "\n")
	}
}

.showOptions <- function(opt) {
	nms <- c("memfrac", "tempdir", "datatype", "progress", "todisk", "verbose") 
	for (n in nms) {
		v <- eval(parse(text=paste0("opt$", n)))
		cat(paste0(substr(paste(n, "         "), 1, 10), ": ", v, "\n"))
	}
}
 
terraOptions <- function(...) {
	dots <- list(...)
	opt <- .terra_environment$options@ptr
	if (length(dots) == 0) {
		.showOptions(opt)
	} else {
		opt <- .setOptions(opt, dots)			
		.terra_environment$options@ptr <- opt
	}
}


tmpFiles <- function(old=FALSE, remove=FALSE) {
	d <- .terra_environment$options@ptr$tempdir
	if (old) {
		f <- list.files(dirname(d), recursive=TRUE, pattern="^spat_", full.names=TRUE)
	} else {
		f <- list.files(d, pattern="^spat", full.names=TRUE)
	}
	if (remove) {
		file.remove(f) 
		return(invisible(f))
	} else {
		return(f)
	}
}


