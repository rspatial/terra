# Author: Robert J. Hijmans
# Date : October 2018
# Version 1.0
# License GPL v3

show_messages <- function(x, f="") {
	if (x@ptr$messages$has_warning) { 
		messages <- paste0(f, ": ", paste(x@ptr$messages$warnings, collapse="\n"))
		x@ptr$messages$warnings <- ""
		x@ptr$messages$has_warning <- FALSE
		warning(messages, call.=FALSE)
	}
	if (x@ptr$messages$has_error) {
		emsg <- x@ptr$messages$error
		x@ptr$messages$error <- ""
		x@ptr$messages$has_error <- FALSE	
		stop(paste0("(", f, ") ", emsg), call.=FALSE)
	}
	x
}

