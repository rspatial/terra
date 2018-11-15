
.terra_environment <- new.env()

 
.create_options <- function() {
	opt <- methods::new("SpatOptions")
	opt@ptr <- SpatOptions$new()
	opt@ptr$tempdir <- tempdir()
	.terra_environment$options <- opt
}
 
.setOptions <- function(x, opt) {
	nms <- names(opt)
	s <- nms %in% c("tempdir", "memfrac", "datatype", "filetype")
	
	if (any(!s)) {
		warning(paste(nms[!s], collapse = ", "), " is/are not valid or cannot be set")
	}
		
	if (any(s)) {
		i <- nms %in% c("datatype", "filetype")
		nms[i] <- paste0("def_", nms[i])
		for (i in which(s)) {
			opt[[nms[i]]] <- opt[[i]]
		}
	}
	opt
} 
 
.runOptions <- function(filename="", overwrite=FALSE, writeopt=list()) {
	writeopt$x = SpatOptions$copy(.terra_environment$options@ptr)
	writeopt$filename = filename
	writeopt$overwrite = overwrite
	opt <- do.call(.setOptions, writeopt)
	#show_messages(opt)
	opt
}


 
terraOptions <- function(...) {

	dots <- list(...)
	opt <- .terra_environment$options@ptr
	
	if (length(dots) == 0) {
		cat("options for package terra\n")
		cat("memfrac     :" , opt$memfrac, "\n")
		cat("tempdir     :" , opt$tempdir, "\n")
		cat("datatype    :" , opt$datatype, "\n")
		cat("filetype    :" , opt$filetype, "\n")
		
	} else {

		opt <- .setOptions(opt, dots)			
		.terra_environment$options@ptr <- opt

	}
}

