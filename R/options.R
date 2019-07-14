
.terra_environment <- new.env()

 
.create_options <- function() {
	opt <- methods::new("SpatOptions")
	opt@ptr <- SpatOptions$new()
	opt@ptr$tempdir <- tempdir()
	.terra_environment$options <- opt
}
 
.setOptions <- function(x, opt) {
	nms <- names(opt)
	s <- nms %in% c("progress", "tempdir", "memfrac", "datatype", "filetype", "filename", "overwrite")
	
	if (any(!s)) {
		warning(paste(nms[!s], collapse = ", "), " invalid option(s)")
	}
		
	if (any(s)) {
		for (i in which(s)) {
			x[[nms[i]]] <- opt[[i]]
		}
	}
	x
} 
 
.runOptions <- function(filename="", overwrite=FALSE, wopt=list()) {
	if (!is.list(wopt)) {
		stop("wopt must be a list")
	}
	ptr <- .terra_environment$options@ptr
	opt <- ptr$copy(ptr)
	
	if (!is.null(unlist(wopt))) {
		wopt$filename <- filename[1]
		wopt$overwrite <- overwrite[1]
		opt <- .setOptions(opt, wopt)
	} else {
		opt$filename <- filename[1]
		opt$overwrite <- overwrite[1]
	}
	#show_messages(opt)
	opt
}

.showOptions <- function(opt) {
	cat("Options for package 'terra'\n")
	cat("memfrac     :" , opt$memfrac, "\n")
	cat("tempdir     :" , opt$tempdir, "\n")
	cat("datatype    :" , opt$def_datatype, "\n")
	cat("filetype    :" , opt$def_filetype, "\n")
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


