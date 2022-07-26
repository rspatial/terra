

.orphanTmpFiles <- function() {

	objects <- ls(envir=globalenv())
	ftmp <- list()
	for (i in seq_along(objects)) {
		x <- get(objects[i], envir=globalenv())
		if (inherits(x, "SpatRaster")) {
			ftmp[[i]] <- sources(x)
		}
	}
	ftmp <- unique(unlist(ftmp))
	ftmp <- ftmp[ftmp != ""]
	pattrn <- "^spat_.*tif$"
	i <- grep(pattrn, basename(ftmp))
	ftmp <- ftmp[i]
	ff <- list.files(tempdir(), pattern=pattrn, full.names=TRUE)
	i <- !(basename(ff) %in% basename(ftmp))
	ff[i]

}




tmpFiles <- function(current=TRUE, orphan=FALSE, old=FALSE, remove=FALSE) {

	if (!(old | current | orphan)) {
		error("tmpFiles", "at least one of 'orphan', 'current' and 'old' must be set to TRUE")
	}

	opt <- spatOptions()
	d <- opt$tempdir
	f <- NULL
	if (old) {
		if (normalizePath(tempdir()) != normalizePath(d)) {
			warn("tmpFiles", "old files can only be found if terra uses the R tempdir")
		} else {
			f <- list.files(dirname(d), recursive=TRUE, pattern="^spat_", full.names=TRUE)
			f <- grep("Rtmp", f, value=TRUE)
			if ((length(f) > 0) && (!current)) {
				i <- grep(d, f)
				if (length(i) > 0) {
					f <- f[-i]
				}
			}
		}
	}

	if (current) {
		ff <- list.files(d, pattern="^spat", full.names=TRUE)
		f <- c(f, ff)
	} else if (orphan) {
		fo <- .orphanTmpFiles()
		f <- c(f, fo) # for if old=TRUE
	}


	if (remove) {
		file.remove(f)
		return(invisible(f))
	} else {
		return(f)
	}
}

