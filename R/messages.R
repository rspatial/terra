# Author: Robert J. Hijmans
# Date : October 2018
# Version 1.0
# Licence GPL v3

.messages <- function(x, f="") {
	if (x@ptr$warning) { 
		messages <- paste0(f, ": ", paste(x@ptr$warning_message, collapse="\n"))
		x@ptr$warning_message <- ""
		x@ptr$warning <- FALSE
		warning(messages, call.=FALSE)
	}
	if (x@ptr$error) {	
		stop(paste0(f, ": ", x@ptr$error_message), call.=FALSE)
	}			
}

