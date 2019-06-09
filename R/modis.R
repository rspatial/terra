# Author: Robert J. Hijmans
# Date : June 2019
# Version 1.0
# License GPL v3


modisqc <- function(x, nbits, se, reject, filename="", overwrite=FALSE, wopt=list(), ... ) {
	opt <- .runOptions(filename, overwrite, wopt)
	x@ptr <- x@ptr$modisqc(nbits, se, reject, opt) 
	show_messages(x)
}

