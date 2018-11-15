
.terra_environment <- new.env()

 
.create_options <- function() {
	opt <- methods::new("SpatOptions")
	opt@ptr <- SpatOptions$new()
	opt@ptr$tempdir <- tempdir()
	.terra_environment$options <- opt
}
 
.runOptions <- function(filename="", filetype="", datatype="", overwrite=FALSE) {
	opt <- SpatOptions$new(.terra_environment$options@ptr)
	opt$filename = filename;
	opt$filetype = filetype;
	opt$datatype = datatype;
	opt$overwrite = overwrite;
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
	
		nms <- names(dots)
		s <- nms %in% c("tempdir", "memfrac", "datatype", "filetype")
		
		if (any(!s)) {
			warning(paste(nms[!s], collapse = ", "), " is/are not valid or cannot be set")
		}
		
		if (any(s)) {
			i <- nms %in% c("datatype", "filetype")
			nms[i] <- paste0("def_", nms[i])
			for (i in which(s)) {
				opt[[nms[i]]] <- dots[[i]]
			}
			.terra_environment$options@ptr <- opt
		}
	}
}

