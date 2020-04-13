# Author: Robert J. Hijmans
# Date : October 2018
# Version 1.0
# License GPL v3

show_messages <- function(x, f="") {
	if (methods::.hasSlot(x, "ptr")) {
		if (x@ptr$messages$has_warning) { 
			messages <- paste0(f, ": ", paste(x@ptr$messages$warnings, collapse="\n"))
			x@ptr$messages$warnings <- ""
			x@ptr$messages$has_warning <- FALSE
			warning(paste(messages, collapse="\n"), call.=FALSE)
		}
		if (x@ptr$messages$has_error) {
			emsg <- x@ptr$messages$error
			x@ptr$messages$error <- ""
			x@ptr$messages$has_error <- FALSE	
			stop(paste0("(", f, ") ", emsg), call.=FALSE)
		}
		return(x)
	} else {
		if (x$messages$has_warning) { 
			messages <- paste0(f, ": ", paste(x$messages$warnings, collapse="\n"))
			x$messages$warnings <- ""
			x$messages$has_warning <- FALSE
			warning(paste(messages, collapse="\n"), call.=FALSE)
		}
		if (x$messages$has_error) {
			emsg <- x$messages$error
			x$messages$error <- ""
			x$messages$has_error <- FALSE	
			stop(paste0("(", f, ") ", emsg), call.=FALSE)
		}
		return(x)
	}
}

