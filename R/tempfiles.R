
orphanTmpFiles <- function() {
	objects <- ls(env=globalenv())
	ftmp <- list()
	for (i in seq_along(objects)) {
		x <- get(objects[i], env=globalenv())
		if (inherits(x, "SpatRaster")) {
			ftmp[[i]] <- sources(x)$source
		}
	}
	ftmp <- unique(unlist(ftmp))
	ftmp <- ftmp[ftmp != ""]
	pattrn <- "^spat_.*tif$"
	i <- grep(pattrn, basename(ftmp))
	ftmp <- ftmp[i]
	ff <- list.files(tempdir(), pattern=pattrn, full=TRUE)
	i <- !(basename(ff) %in% basename(ftmp))
	ff[i]
}

