
staleTmpFiles <- function() {
	pattrn <- "^spat_.*tif$"
	ff <- list.files(dirname(tempdir()), pattern="^spat_.*\\.tif$", full.names=TRUE, recursive=TRUE)
	grep("Rtmp", ff, value=T)
}


orphanTmpFiles <- function() {

	objects <- ls(envir=globalenv())
	ftmp <- list()
	for (i in seq_along(objects)) {
		x <- get(objects[i], envir=globalenv())
		if (inherits(x, "SpatRaster")) {
			ftmp[[i]] <- sources(x)$source
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

