# Author: Robert J. Hijmans
# Date : October 2018
# Version 1.0
# License GPL v3

error <- function(f, emsg="", ...) {
	stop("[", f, "] ", emsg, ..., call.=FALSE)
}

warn <- function(f, wmsg="", ...) {
	warning("[", f, "] ", wmsg, ..., call.=FALSE)
}

messages <- function(x, f="") {
	#g <- gc(verbose=FALSE)
	if (methods::.hasSlot(x, "ptr")) {
		if (x@ptr$has_warning()) {
			warn(f, paste(unique(x@ptr$getWarnings()), collapse="\n"))
		}
		if (x@ptr$has_error()) {
			error(f, x@ptr$getError())
		}
	} else {
		if (x$has_warning()) {
			warn(f, paste(unique(x$getWarnings()), collapse="\n"))
		}
		if (x$has_error()) {
			error(f, x$getError())
		}
	}
	x
}


mem_info <- function(x, n=1) {
	#print=TRUE
	n <- max(0,n)
	opt <- spatOptions()
	opt$ncopies = n;
	v <- x@ptr$mem_needs(opt)
	memmin <- opt$memmin
	memmax <- opt$memmax
	#if (print) {
		gb <- 1024^3 / 8  #
		cat("\n------------------------")
		cat("\nMemory (GB) ")
		cat("\n------------------------")

		cat(paste("\ncheck threshold :", opt$memmin / gb, "(memmin)"))
		if (memmax > 0) {
			cat(paste("\navailable       :",  round(v[2] / gb, 2), "(memmax)"))
		} else {
			cat(paste("\navailable       :",  round(v[2] / gb, 2)))
		}
		cat(paste0("\nallowed (", round(100* v[3]) , "%)   : ", round(v[3] * v[2] / gb, 2)))

		cat(paste0("\nneeded (n=", n, ")   ", ifelse(n<10, " : ", ": "), round(v[1] / gb, 2)))
		cat("\n------------------------")
		cat(paste("\nproc in memory  :", round(v[5]) != 0))
		cat(paste("\nnr chunks       :", ceiling(nrow(x)/v[4])))
		cat("\n------------------------\n")
	#}
	names(v) <- c("needed", "available", "memfrac", "chunksize")
	invisible(v)
}


free_RAM <- function() {
	opt <- spatOptions()
	x <- rast()
	v <- x@ptr$mem_needs(opt)
	v[2] / 128
}



