# Author: Robert J. Hijmans
# Date : October 2018
# Version 1.0
# License GPL v3

show_messages <- function(x, f="") {
	if (methods::.hasSlot(x, "ptr")) {
		if (x@ptr$messages$has_warning) { 
			warns <- x@ptr$messages$getWarnings()
			warning(paste(warns, collapse="\n"), call.=FALSE)
		}
		if (x@ptr$messages$has_error) {
			emsg <- x@ptr$messages$getError()
			stop(paste0("[", f, "] ", emsg), call.=FALSE)
		}
		return(x)
	} else {
		if (x$messages$has_warning) { 
			messages <- paste0(f, ": ", paste(x$messages$getWarnings(), collapse="\n"))
			warning(paste(messages, collapse="\n"), call.=FALSE)
		}
		if (x$messages$has_error) {
			emsg <- x$messages$getError()
			stop(paste0("[", f, "] ", emsg), call.=FALSE)
		}
		return(x)
	}
}


.mem_info <- function(x, n=1, print=TRUE) {
	n <- max(0,n)
	opt <- .getOptions()	
	v <- x@ptr$mem_needs(n, opt)
	if (print) {
		gb <- 1073741824 
		cat("\n----------------------")
		cat("\nMemory (GB) ")
		cat("\n----------------------")
		cat(paste("\navailable     :", round(v[2] / gb, 2)))
		cat(paste0("\nallowed (", round(100* v[3]) , "%) : ", round(v[3] * v[2] / gb, 2)))
		cat(paste("\nneeded        :", round(v[1] / gb, 2)))
		cat("\n----------------------")
		cat(paste("\nproc in memory:", round(v[5]) != 0))
		cat(paste("\nnr chunks     :", ceiling(nrow(x)/v[4])))
		cat("\n----------------------\n")
	} else {
		names(v) <- c("needed", "available", "memfrac", "chunksize")
		v
	}
}
	
