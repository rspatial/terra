
.terra_environment <- new.env()

 
.create_options <- function() {
	opt <- methods::new("SpatOptions")
	opt@ptr <- SpatOptions$new()
	opt@ptr$tempdir <- tempdir()
	.terra_environment$options <- opt
}
 
 
terraOptions <- function(...) {

	dots <- list(...)
	opt <- .terra_environment$options@ptr
	
	if (length(dots) == 0) {
		cat("options for package terra\n")
		cat("tempdir     :" , opt$tempdir, "\n")
		cat("memfrac     :" , opt$memfrac, "\n")
		cat("datatype    :" , opt$datatype, "\n")
		
	} else {
	
		nms <- names(dots)
		s <- nms %in% c("tempdir", "memfrac", "datatype")
		
		if (any(!s)) {
			warning(paste(nms[!s], collapse = ", "), " is/are not valid or cannot be set")
		}
		
		if (any(s)) {
			for (i in which(s)) {
				opt[[nms[i]]] <- dots[[i]]
			}
			.terra_environment$options@ptr <- opt
		}
	}
}

